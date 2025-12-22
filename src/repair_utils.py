import os

def run_getTargetArray(ref_file,cen_list_file,outdir,outdir_cen_array):
    ref_asm = ref_file 
    with open(cen_list_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            array_name = items[0] + ':' + items[1] + '-' + items[2]
            cmd = 'samtools faidx ' + ref_asm + ' ' + array_name + ' > ' +  outdir_cen_array + '/' + array_name + '.tmp.fa'
            os.system(cmd)
            cmd = 'seqkit seq -w 0 ' + outdir_cen_array + '/' + array_name + '.tmp.fa' + ' > ' + outdir_cen_array + '/' + array_name + '.fa'
            os.system(cmd)
            cmd = 'rm ' +  outdir_cen_array + '/' + array_name + '.tmp.fa'
            os.system(cmd)
    cmd = 'cat ' + outdir_cen_array + '/*.fa > ' + outdir + '/all.cenarray.fa'
    os.system(cmd)


def getFilterRecord(paf_file,cen_list_file,outfile,MAPQ_thr,break_contig_merge_thr,min_array_thr,overlap_merge_thr): 
    target_arrays = {}
    
    # 构建最初的target_arrays
    with open(cen_list_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            array_name = items[0] + ':' + items[1] + '-' + items[2]
            target_arrays[array_name] = {}
    
    # 从paf 文件中读取，保存到target_arrays
    # key1: array, key2:hifiasm-contig_strand, value: 对应记录
    with open(paf_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            # print(items)
            if int(items[11]) <= MAPQ_thr:
                continue 
            if items[0] in target_arrays.keys():
                if items[5] + '_' +  items[4] not in target_arrays[items[0]].keys():
                    target_arrays[items[0]][items[5] + '_' +  items[4] ] = []
                    target_arrays[items[0]][items[5] + '_' +  items[4] ].append([items[0],items[1],int(items[2]),int(items[3]),'+',items[5],int(items[7]),int(items[8]), items[4],items[11]])
                else:
                    target_arrays[items[0]][items[5] + '_' +  items[4] ].append([items[0],items[1],int(items[2]),int(items[3]),'+',items[5],int(items[7]),int(items[8]), items[4],items[11]])
    
    
    
    # 对于记录按照hifiasm start进行排序
    sorted_target_arrays = {}
    for i in target_arrays.keys():
        sorted_target_arrays[i] = {}
        for j in target_arrays[i].keys():
            sorted_target_arrays[i][j] = sorted(target_arrays[i][j],key=lambda x:x[6])

    # 对于每一个sorted_target_arrays[key1][key2]，合并其中记录基于，大小break_contig_merge_thr
    target_arrays_tmp = {}
    for i in sorted_target_arrays.keys():
        target_arrays_tmp[i] = {}
        for j in sorted_target_arrays[i].keys():
            target_arrays_tmp[i][j] = []
            if len(sorted_target_arrays[i][j]) == 1:
                target_arrays_tmp[i][j] = sorted_target_arrays[i][j]
            else: 
                merge_records = []
                init_record = sorted_target_arrays[i][j][0]
                for k in range(len(sorted_target_arrays[i][j]) -1):
                    next_record = sorted_target_arrays[i][j][k + 1]
                    if abs(next_record[6] - init_record[7]) < break_contig_merge_thr:
                        if next_record[7] > init_record[7]:
                            # 进行合并
                            init_record[7] = next_record[7]
                            if init_record[2] > next_record[2]:
                                init_record[2] = next_record[2]
                            if init_record[3] < next_record[3]:
                                init_record[3] = next_record[3]
                        else:
                            continue
                    else:
                        merge_records.append(init_record)
                        init_record = next_record
                merge_records.append(init_record)

                target_arrays_tmp[i][j] = merge_records
  

    target_arrays = {}
    # 重新生成记录，并去掉小于100k的
    # 重新存成key: array, value: 记录
    for i in target_arrays_tmp.keys():
        target_arrays[i] = []
        for j in target_arrays_tmp[i].keys():
            for k in target_arrays_tmp[i][j]:
                if k[7] - k[6] > min_array_thr:
                    target_arrays[i].append(k)

    # second contig必须唯一对应target array,如果不唯一，尝试选择比对长度最长的 
    # note: *着丝粒可行，contig唯一对应，但是泛化到其他区间，例如SD可能不行, 已修改
    outfile = open(outfile,'w')
    filter_target_array = {}
    
    
    for i in target_arrays.keys():
        # if i != 'chr10_RagTag_hap1:114448959-115450181':
        #     continue
        sorted_regions = sorted(target_arrays[i],key=lambda x:x[2]) # 按照target asm array的start排序
        # 去重，后一个如果跟前一个重叠，且重叠大于后者的try 30%，去掉后者
        final_pairs = []
        if len(sorted_regions) < 2:
            final_pairs = sorted_regions
            if len(sorted_regions) == 0:
                print(i)
        else:
            final_pairs.append(sorted_regions[0])
            for j in range(len(sorted_regions) - 1):
                if sorted_regions[j + 1][2] > final_pairs[-1][3]:
                    final_pairs.append(sorted_regions[j + 1])
                else:
                    if sorted_regions[j + 1][3] <= final_pairs[-1][3]:
                        continue
                    else:
                        # 检查重叠大小
                        if (final_pairs[-1][3]- sorted_regions[j + 1][2]) / int(sorted_regions[j + 1][1]) > overlap_merge_thr:
                            continue
                        else:
                            final_pairs.append(sorted_regions[j + 1])
        
        # final_pairs, 输出了最终的记录
        contig_table = {} # 建立了一个key:hifiasm contig, value:记录
        for j in final_pairs:
            if j[5] not in contig_table.keys():
                contig_table[j[5]] = [j]
            else:
                contig_table[j[5]].append(j)
  
        
        filter_target_array[i] = {}
        # 遍历contig_table，如果有多个记录保留一个最大的一个，如果该hifiasm contig在该array上是碎裂的，就选择最大的主体
        for j in contig_table.keys():
            record = contig_table[j] # [verkko cen_array,array lenth,array_start,array_end,'+',hifiasm contig name,start,end, hifiasm strand, MQ]
            if len(record) == 1:
                filter_target_array[i][j] = record[0]
            else:
                # 找最长的
                max_length = -1
                max_record = []
                for k in record:
                    if (int(k[7]) - int(k[6])) > max_length:
                        max_length = (int(k[7]) - int(k[6]))
                        max_record = k
                filter_target_array[i][j] = max_record
    
    # 进行过滤操作
    # final_match key:hifiasm-contig, value:array
    final_target_array = {} 
    # 遍历filter_target_array，确定记录的hifiasm-contig是verkko的array之后才加入
    for i in filter_target_array.keys():
        final_target_array[i] = []    
        for j in filter_target_array[i].keys():
            final_target_array[i].append(filter_target_array[i][j])
    
    for i in final_target_array.keys():
        for j in final_target_array[i]:
            record = j
            outfile.write(record[0] + '\t' + record[1] + '\t' + str(record[2]) + '\t' + str(record[3]) + '\t' + record[4] + '\t' + record[5] + '\t' + str(record[6]) + '\t' + str(record[7]) + '\t' +record[8] +'\t' + record[9]+ '\n')
    outfile.close()
 

  
def buildErrorRegions(regions_file,error_dir): 
    if not os.path.exists(error_dir):
        os.mkdir(error_dir)
    arrays = {}
    with open(regions_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            start = int(items[0].split(':')[1].split('-')[0])
            array = items[0]
            if items[-1] == 'HET':
                continue
            if array not in arrays:
                    arrays[array] = []
                    arrays[array].append([int(items[1]) - start, int(items[2]) - start])
            else:
                if items[-1] == 'HET':
                    continue
                arrays[array].append([int(items[1]) - start, int(items[2]) - start])
            
    for i in arrays.keys():
        outfile = error_dir + '/' + i + '.xls'
        outfile = open(outfile,'w')
        for j in arrays[i]:
            outfile.write(str(j[0]) + '\t' + str(j[1]) + '\n')
        outfile.close()

def buildInitNucflagRegions(matching_records_file,seconda_asm_array_bed,outrepair_table):
    outseconda_asm_array_bed = open(seconda_asm_array_bed,'w')
    outrepair_table = open(outrepair_table,'w')
    
    repair_record = {}
    
    with open(matching_records_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            array = items[0]
            # 不能用contig name作为key
            arrayinfo = array.split(':')
            region = arrayinfo[1].split('-')
            contig = items[5]
            start = int(items[6]) + 1
            end = int(items[7])
            outseconda_asm_array_bed.write(contig + '\t' + str(start) + '\t' + str(end) + '\t' + str(end - start) +  '\n')
            if array not in repair_record.keys():
                repair_record[array] = [[region[0],region[1],'+',contig,str(start),str(end),items[8]]]
            else:
                repair_record[array].append([region[0],region[1],'+',contig,str(start),str(end),items[8]])
                
    outseconda_asm_array_bed.close()
    
    for i in repair_record.keys():
        count = 1
        contig = i.split(':')[0]
        for j in repair_record[i]:
            outrepair_table.write(contig + '\t' + j[0] + '\t' + j[1] + '\t' + j[2] + '\t' + j[3] + '\t' + j[4] + '\t' + j[5] + '\t' + j[6] + '\t' + str(count) + '\n')
            count += 1
    
    outrepair_table.close()
    
def run_nucflag(hifi_reads_files,target_ref_file,out_nucflag_dir,hifi_reads_mapping_threads,nucflag_ignore_regions,nucflag_config,coordinates_file,error_dir):
    # mapping hifi reads

    bam_files_list = []
    
    if not os.path.exists(out_nucflag_dir + '/nucflag.bam.bai'): # if bam already there, skip alignment
        inter_bam_dir = out_nucflag_dir + '/temp_bam'
        if not os.path.exists(inter_bam_dir):
            os.mkdir(inter_bam_dir)
    
        for i in hifi_reads_files:
            # mapping
            # pbmm2 align         --log-level DEBUG         --preset SUBREAD         --min-length 5000         -j 24 ref reads
            filename = i.split('/')[-1]
            cmd = 'pbmm2 align --log-level DEBUG --preset SUBREAD --min-length 5000 -j ' + \
                str(hifi_reads_mapping_threads) + ' ' + target_ref_file + ' ' + i + ' > ' + inter_bam_dir + '/' + filename + '.init.bam'
            os.system(cmd)
            # samtools view -F 2308 -@  -u 
            cmd = 'samtools view -F 2308 -@ ' + str(hifi_reads_mapping_threads) + ' -u ' + inter_bam_dir + '/' + filename + '.init.bam' + \
                ' -o ' + inter_bam_dir + '/' + filename +'.filtered.bam'
            os.system(cmd)
            # samtools sort -T temp -m 4G -@ 24 -
            cmd = 'samtools sort -T temp -m 4G -@ ' +str(hifi_reads_mapping_threads) + \
                ' -o ' + inter_bam_dir + '/' + filename + '.bam' +' ' + inter_bam_dir + '/' + filename +'.filtered.bam'
            os.system(cmd)
            
            bam_files_list.append(inter_bam_dir + '/' + filename + '.bam')
            
            # rm init and after-filter bam, keep after-sort bam
            cmd = 'rm ' + inter_bam_dir + '/' + filename + '.init.bam'
            os.system(cmd)
            cmd = 'rm ' + inter_bam_dir + '/' + filename + '.filtered.bam'
            os.system(cmd)
            
        # merge and remove    
        bam_files = ''
        for i in bam_files_list:
            bam_files += i + '  '
    
        # samtools merge -@ 24 - bam1 bam2 bam3
        cmd = 'samtools merge -@ ' + str(hifi_reads_mapping_threads) + ' ' + ' -o ' + out_nucflag_dir + '/nucflag.init.bam'  + '  ' + bam_files
        os.system(cmd)
        # samtools sort -T temp -m 4G -@ 24 -;} > results/nucflag/NA18534_hifi.bam
        cmd = 'samtools sort -T temp -m 4G -@ ' + str(hifi_reads_mapping_threads) + ' ' + out_nucflag_dir + '/nucflag.init.bam' + ' > ' +  out_nucflag_dir + '/nucflag.bam' 
        os.system(cmd)
        # rm
        cmd = 'rm ' + out_nucflag_dir + '/nucflag.init.bam'
        os.system(cmd)
        cmd = 'rm -rf ' + inter_bam_dir
        os.system(cmd)
        
        # index 
        cmd = 'samtools index ' + out_nucflag_dir + '/nucflag.bam' 
        os.system(cmd)
    
    # run nucflag
    cmd = 'nucflag  -i ' + out_nucflag_dir + '/nucflag.bam'  + \
        ' -d ' + out_nucflag_dir + '/png' + '  ' + \
            ' -o ' + out_nucflag_dir + '/nucflag.misassemblies.bed' + ' ' + \
                ' -t ' + str(hifi_reads_mapping_threads) + ' ' + \
                ' -p ' + str(hifi_reads_mapping_threads) + ' ' + ' -s ' + out_nucflag_dir + '/nucflag.status.bed' + ' ' + \
                    ' -c ' +  nucflag_config + ' ' 
    if coordinates_file != '':
        cmd += ' -b ' +  coordinates_file + ' ' 
    
    if nucflag_ignore_regions != '':
        cmd += ' --ignore_regions ' + nucflag_ignore_regions
        
    
    os.system(cmd)
    
    # get error regions
    if coordinates_file != '':
        buildErrorRegions(out_nucflag_dir + '/nucflag.misassemblies.bed',error_dir)
    
    

def run_repair(script, match_file,target_cen_array,second_asm_array,error_dir,second_error_dir,outrepairdir,kmer_size=5000,kmer_number=2000,extend_length=1000):

    contigs_repair_record = {}
    with open(match_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            name = items[0] + ':' + items[1] + '-' + items[2] + '@' + items[3]
            if name not in contigs_repair_record.keys():
                contigs_repair_record[name] = [line]
            else:
                contigs_repair_record[name].append(line)
    
    for name in contigs_repair_record.keys():
        out_repair_file_name = outrepairdir + '/' + name + '.repair.xls'
        out_repair_file = open(out_repair_file_name,'w')
        
        for j in contigs_repair_record[name]:
            out_repair_file.write(j + '\n')
        out_repair_file.close()
        
    
    for name in contigs_repair_record.keys():
        out_repair_file_name = outrepairdir + '/' + name + '.repair.xls'
        region = name.split('@')[0]
        error_regions_file = error_dir + '/' + region + '.xls'
        if not os.path.exists(error_regions_file):
            error_regions_file = '-1'
        error_regions_flag = name.split('@')[1]
        
        target_fa_file = target_cen_array + '/' + region + '.fa'
        
        outfa_file = outrepairdir + '/' + name + '.repaired.fa'
        outlog = outrepairdir + '/' + name + '.repaired.log'
        print(outlog)
        cmd = 'python ' + script + ' -r ' + out_repair_file_name + \
            ' -ts ' + target_fa_file + ' -tf ' + error_regions_flag + ' -te ' + error_regions_file + \
                ' -ssd ' + second_asm_array + ' -sed ' + second_error_dir + \
                    ' -omd ' + outrepairdir + ' -os ' + outfa_file + ' -ol ' + outlog + \
                        ' -k ' + str(kmer_size) + ' -kn ' + str(kmer_number) + ' -et ' + str(extend_length)
        
        os.system(cmd)

def getContigLength(fai_file):
    contig_length = {}
        
    with open(fai_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            contig_length[items[0]] = int(items[1])

    return contig_length

def buildAllErrorRegions(error_file,fai_file,outfile,merge_distance):
    contig_length = getContigLength(fai_file)
    chr_bed = {}
    with open(error_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            if 'HET' in line:
                continue
            chr = items[0]
            start = int(items[1])
            end = int(items[2])
            if chr not in chr_bed.keys():
                chr_bed[chr] = []
                chr_bed[chr].append([start,end])
            else:
                chr_bed[chr].append([start,end])
                
    chr_merge_bed = {}
    for i in chr_bed.keys():
        chr_merge_bed[i] = []
        if len(chr_bed[i]) == 1:
            chr_merge_bed[i]  = chr_bed[i]
        else:
            sorted_bed = sorted(chr_bed[i],key=lambda x:x[0])
            # 先去掉重复,如果下一个的end小于上一个，则删掉
            remove_dup = []
            remove_dup.append(sorted_bed[0])
            for j in range(len(sorted_bed) - 1):
                if sorted_bed[j + 1][1] <= sorted_bed[j][1]:
                    continue
                else:
                    remove_dup.append(sorted_bed[j + 1])
            # 合并
            init_region = remove_dup[0]
            for j in range(len(remove_dup) - 1):
                if (remove_dup[j + 1][0] - init_region[1]) < merge_distance:
                    init_region = [init_region[0], remove_dup[j+1][1]]
                else:
                    chr_merge_bed[i].append(init_region)
                    init_region = remove_dup[j + 1]
            chr_merge_bed[i].append(init_region)
    final_regions = []
    
    for i in chr_merge_bed.keys():
        for j in chr_merge_bed[i]:
            start = int(j[0]) - int(merge_distance / 2)
            if start < 1:
                start = 1
            end = int(j[1]) + int(merge_distance / 2)
            if end > contig_length[i]:
                end = contig_length[i]
            
            final_regions.append([i, start, end])
    
    outfile = open(outfile,'w')
    for i in final_regions:
        outfile.write(i[0] + '\t' + str(i[1]) + '\t' + str(i[2]) + '\n')
    outfile.close()
        

