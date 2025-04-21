
import random
import argparse
import os

def readErrorRegions(file,flag,array_length):
    # 在这里进行错误区间的读取
    # 格式start end (坐标应为相对坐标)
    error_regions = []
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            error_regions.append([int(items[0]),int(items[1])])
    error_regions = sorted(error_regions,key=lambda x:x[0])
    
    if flag == '-':
        new_error_regions = []
        for j in error_regions[::-1]:
            new_error_regions.append([array_length - 1 - j[1],array_length - 1 - j[0]])
        error_regions = new_error_regions
    
    return error_regions    

def get_reverse_complement(dna_sequence):
    # 定义碱基互补对
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
    # 生成互补链
    complement_sequence = ''.join(complement[base] for base in dna_sequence)
    # 反转互补链得到反链
    reverse_complement_sequence = complement_sequence[::-1]
    return reverse_complement_sequence

def readSeq(file,flag):
    # 要求输入序列为一条，不要format        
    seq = ''
    with open(file,'r') as f:
        line = f.readline()[:-1]
        seq = f.readline()[:-1]
    if flag == '-':
        seq = get_reverse_complement(seq)    
    
    return seq
            

def exchange_regions(error_regions,array_length,extend_length):
    merge_exchange_regions = [] # 保存需要交换的区间
    if len(error_regions) == 1: # 当只有一个时候，不需要合并，直接延展即可
        exchange_start = error_regions[0][0] - extend_length 
        if exchange_start < 0: # 注意不要超过array的start
            exchange_start = 0 # 超过了就设置为array start
        exchange_end = error_regions[0][1] + extend_length
        if exchange_end > array_length - 1: #  注意不要超过array的end
            exchange_end = array_length - 1 # 超过了就设置为array end
        merge_exchange_regions.append([exchange_start,exchange_end])
    else: # 超过1个
        init_exchange_region = error_regions[0] # 初始
        for i in range(len(error_regions) - 1):
            # ----()()-----
            if error_regions[i+1][0] - init_exchange_region[1] <= extend_length: # 当后一个和前一个的距离小于等于
                init_exchange_region[1] = error_regions[i+1][1] # 进行区间合并
            else: # 否则进行输出
                exchange_start = init_exchange_region[0] - extend_length 
                if exchange_start < 0:
                    exchange_start = 0
                exchange_end = init_exchange_region[1] + extend_length
                if exchange_end > array_length - 1:
                    exchange_end = array_length - 1
                merge_exchange_regions.append([exchange_start,exchange_end])
                init_exchange_region = error_regions[i+1]
                
        exchange_start = init_exchange_region[0] - extend_length 
        if exchange_start  < 0:
            exchange_start = 0
        exchange_end = init_exchange_region[1] + extend_length
        if exchange_end > array_length - 1:
            exchange_end = array_length - 1
        merge_exchange_regions.append([exchange_start,exchange_end])
    return merge_exchange_regions
    
def buildKmeMatch(targetasm_fa,secondasm_fa, kmer_size,kmer_number):
    all_kmers = {}
    for i in range(len(targetasm_fa) - kmer_size):
        kmer = targetasm_fa[i:i+kmer_size]
        if kmer not in all_kmers.keys():
            all_kmers[kmer] = [i]
        else:
            all_kmers[kmer].append(i)
    # 去除>=2的kmer
    uniq = 0
    multiple = 0
    uniq_kmer_in_targetasm = {}
    for i in all_kmers.keys():
        if len(all_kmers[i]) >= 2:
            multiple += 1
        else:
            uniq += 1
            uniq_kmer_in_targetasm[i] = all_kmers[i]

    # 随机筛选kmer size个
    random_kmers = random.sample(list(uniq_kmer_in_targetasm.keys()), kmer_number)
    random_kmers_set = set()
    for i in random_kmers:
        random_kmers_set.add(i)
    
    # 然后在secondasm组装中去掉多比对的kmer
    all_kmers_in_secondasm = {} 
    for i in range(len(secondasm_fa) - kmer_size):
        kmer = secondasm_fa[i:i+kmer_size]
        if kmer in random_kmers_set:
            if kmer not in all_kmers_in_secondasm.keys():
                all_kmers_in_secondasm[kmer] = [i]
            else:
                all_kmers_in_secondasm[kmer].append(i)

    # 多比对kmer
    multiple_secondasm_kmers = set()
    for i in all_kmers_in_secondasm.keys():
        if len(all_kmers_in_secondasm[i]) >= 2:
            multiple += 1
            multiple_secondasm_kmers.add(i)
        else:
            uniq += 1

    # 建立最终match列表
    final_match_table = {}
    for i in random_kmers:
        if i in multiple_secondasm_kmers:
            continue
        targetasm_index = uniq_kmer_in_targetasm[i][0]
        secondasm_index = -1
        if i in all_kmers_in_secondasm.keys():
            secondasm_index = all_kmers_in_secondasm[i][0]

        final_match_table[i] = [targetasm_index,secondasm_index]
    # sort
    sorted_final_match_table = sorted(final_match_table.items(),key=lambda x:x[1][0])

    return sorted_final_match_table

def readKmerPair(match_file):
    # 生成block，block内顺序一致
    # 先合并差别10以内的block
    # 放弃最终小于10的block
    # 最后根据顺序进行过滤
    
    sorted_final_match_table = {}
    
    pairs = []

    with open(match_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            pairs.append([int(items[0]),int(items[1]),items[2]])
    blocks = []
    block = []
    for i in range(len(pairs)):
        if pairs[i][1] == -1:
            if len(block) != 0:
                # blocks: [block , block]
                blocks.append(block)
                block = []
        else:    
            # block: [[pair,index],[pair,index],[pair,index]]
            block.append([pairs[i],i])

    if len(block) != 0:
        blocks.append(block)
    
    merge_blocks = []
    if len(blocks) < 2:
        merge_blocks = blocks
    else:
        merge_block = blocks[0]
        for i in range(len(blocks) - 1):
            if blocks[i + 1][0][-1] - blocks[i][-1][-1] < 10:
                merge_block = merge_block + blocks[i+1]
            else:
                merge_blocks.append(merge_block)
                merge_block = blocks[i + 1]
        merge_blocks.append(merge_block)
    
    
    filter_merge_blocks = []
    for i in merge_blocks:
        if len(i) > 10:
            filter_merge_blocks.append(i)
    
    # 再每个block中基于顺序进行pair过滤
    for i in filter_merge_blocks:
        block = i
        pairs = []
        for j in block:
            pairs.append(j[0])
        sorted_final_match_table[pairs[0][2]] = [int(pairs[0][0]),int(pairs[0][1])]
        last = pairs[0][1]
        for j in range(len(pairs) - 1):
            if pairs[j+1][1] - last >= 0:
                sorted_final_match_table[pairs[j+1][2]] = [int(pairs[j+1][0]),int(pairs[j+1][1])]
                last = pairs[j + 1][1]  
    
    return sorted_final_match_table

def fixArray_dev(error_regions,contigs_match_table,targetasm_fa,targetasm_array_length,targetasm_flag):
    exchange_regions = []
    error_regions_repair_flag = {}
    for i in error_regions:
        error_regions_repair_flag[str(i[0]) + '_' + str(i[1])] = 0
    print('error regions:')
    print(error_regions)
    boundary_flag = 0
    repair_flag = 0
    # [updated_sorted_final_match_table,secondasm_fa,secondasm_array_length,secondasm_error_regions]
    for record in contigs_match_table.keys():
        updated_sorted_final_match_table = contigs_match_table[record][0]
       
        for i in error_regions:
            if error_regions_repair_flag[str(i[0]) + '_' + str(i[1])] != 0:
                continue
            
            can_fix = 1
            
            if targetasm_flag == '-': # 当 flag = '-' 时，nucflag结果也得变反
                start_tmp = i[0]
                end_tmp = i[1]
                start = targetasm_array_length - 1 - end_tmp
                end = targetasm_array_length - 1 - start_tmp
            else:
                start = i[0]
                end = i[1]

            targetasm_kmer_index_start = -1
            secondasm_kmer_index_start = -1
            targetasm_kmer_index_end = -1
            secondasm_kmer_index_end = -1
            
            find_start = 0
            find_end = 0
            key_list = list(updated_sorted_final_match_table.keys())
            for j in range(len(key_list)):
                # 遍历找到第一个index超过target start的，则选取上一个pair的index
                # 遍历找到第一个index超过target end的就选取这个pair的index
                
                if find_start != 1:
                    if updated_sorted_final_match_table[key_list[j]][0] > start:
                        if j - 1 < 0:
                            can_fix = 0
                            boundary_flag = 1
                        else:
                            targetasm_kmer_index_start = updated_sorted_final_match_table[key_list[j - 1]][0]
                            secondasm_kmer_index_start = updated_sorted_final_match_table[key_list[j - 1]][1]
                        find_start = 1
                
                if find_end != 1:
                    if updated_sorted_final_match_table[key_list[j]][0] > end:
                        targetasm_kmer_index_end = updated_sorted_final_match_table[key_list[j]][0]
                        secondasm_kmer_index_end = updated_sorted_final_match_table[key_list[j]][1]
                        find_end = 1
                
            if targetasm_kmer_index_start == -1 or targetasm_kmer_index_end == -1:
                can_fix = 0
                boundary_flag = 1
            
            if find_end == 0:
                can_fix = 0
                boundary_flag = 1
            
            if can_fix == 1:
                exchange_regions.append([(targetasm_kmer_index_start,targetasm_kmer_index_end + 1),(record,secondasm_kmer_index_start,secondasm_kmer_index_end + 1)]) # [2,4) [5,7)    
                error_regions_repair_flag[str(i[0]) + '_' + str(i[1])] = 1
                
    print('initial exchange regions:')
    print(exchange_regions)
    final_exchange_regions = []
    for i in exchange_regions:
        overlap_flag = 0
        secondasm_region = i[1] # (secondasm_kmer_index_start,secondasm_kmer_index_end + 1)
        record = secondasm_region[0]
        secondasm_region_start = secondasm_region[1]
        secondasm_region_end = secondasm_region[2]
        secondasm_error_regions = contigs_match_table[record][3]
        
        # 遍历检查区间
        for j in secondasm_error_regions:
            #            ----------
            #      -----             ------------
            if (j[1] <  secondasm_region_start) or (j[0] > secondasm_region_end):
                continue
            else:
                overlap_flag = 1
        if overlap_flag == 0:
            final_exchange_regions.append(i)
    
    print('final exchange regions:')
    print(final_exchange_regions) # Bug: 最终区间需要进行合并，可能会存在重叠性
    
    
    if boundary_flag == 1:
        repair_flag = 2
        
    
    if len(final_exchange_regions) == 0:
        return targetasm_fa,len(targetasm_fa),repair_flag
    
    # 检测，如果targetasm的区间存在重叠，且在两个contig上，则应该将第二个移除，无法修复
    # final_exchange_regions = [[(4475822, 4510359), (('h1tg000014l', '1', '4625766', '-'), 4492198, 4510359)], 
    #                           [(4510400, 4610400), (('h1tg000014l', '1', '4625766', '-'), 4510359, 4710359)], 
    #                           [(4710400, 4810400), (('h1tg000014l', '1', '4625766', '-'), 4892198, 4910359)], 
    #                           [(8540620, 8540630), (('h1tg000041l', '38376304', '43728908', '-'), 3857136, 3878993)]]
    if len(final_exchange_regions) != 1:
        final_exchange_regions_tmp = []
        init_targetasm_region = final_exchange_regions[0][0]
        final_exchange_regions_tmp.append(final_exchange_regions[0])
        for i in range(len(final_exchange_regions) - 1):
            next_targetasm_region = final_exchange_regions[i + 1][0]    
            if next_targetasm_region[0] > final_exchange_regions_tmp[-1][0][1]:
                final_exchange_regions_tmp.append(final_exchange_regions[i + 1])
            else:
                # 如果vekko 的区间重叠，检查重叠是否发生在两个secondasm contigs上，如果是，则删除后面一个
                if final_exchange_regions_tmp[-1][1][0][0] != final_exchange_regions[i + 1][1][0][0]:
                    continue
                else:
                    final_exchange_regions_tmp.append(final_exchange_regions[i + 1])
        print('cross contigs filtering')
        print(final_exchange_regions_tmp)
        final_exchange_regions = final_exchange_regions_tmp
    
    if len(final_exchange_regions) == len(error_regions):
        print('Can fixed')
        repair_flag = 1
    else:
        repair_flag = 3
    
    # Bug: 最终区间需要进行合并，可能会存在重叠性，结合targetasm以及secondasm进行合并
    
    merged_final_exchange_regions = []
    if len(final_exchange_regions) == 1:
        merged_final_exchange_regions = final_exchange_regions
    else:
        init_targetasm_region = final_exchange_regions[0][0]
        init_secondasm_region = final_exchange_regions[0][1]
        for i in range(len(final_exchange_regions) - 1):
            # [(4475822, 8540631), (('h1tg000014l', '1', '4625766', '-'), 4492198, 4510359)]
            next_targetasm_region = final_exchange_regions[i + 1][0] # (4475822, 8540631)
            next_secondasm_region = final_exchange_regions[i + 1][1] # (('h1tg000014l', '1', '4625766', '-'), 4492198, 4510359)
            
            # 目前已经不存在跨contigs的重叠的区间，如果存在重叠一定在同一个contigs上，如果vekko或者secondasm重叠那么更新，将区间合并
            if (next_targetasm_region[0] <= init_targetasm_region[1]) :
                if (next_targetasm_region[1] >  init_targetasm_region[1]):
                    init_secondasm_region = (init_secondasm_region[0],init_secondasm_region[1],next_secondasm_region[2])
                    init_targetasm_region = (init_targetasm_region[0],next_targetasm_region[1])
                else:
                    continue
            elif (next_secondasm_region[1] <= init_secondasm_region[2]) and init_secondasm_region[0][0] == next_secondasm_region[0][0] :
                if (next_secondasm_region[2] >  init_secondasm_region[2]):
                    init_secondasm_region = (init_secondasm_region[0],init_secondasm_region[1],next_secondasm_region[2])
                    init_targetasm_region = (init_targetasm_region[0],next_targetasm_region[1])
                else:
                    continue
            else:
                merged_final_exchange_regions.append([init_targetasm_region,init_secondasm_region])
                init_targetasm_region = final_exchange_regions[i+1][0]
                init_secondasm_region = final_exchange_regions[i+1][1]
        merged_final_exchange_regions.append([init_targetasm_region,init_secondasm_region])
    
    print('merged final exchange regions:')
    print(merged_final_exchange_regions) 
    # [[(4475822, 4610400), (('h1tg000014l', '1', '4625766', '-'), 4492198, 4710359)], 
    #  [(4710400, 4810400), (('h1tg000014l', '1', '4625766', '-'), 4892198, 4910359)], 
    #  [(8540620, 8540630), (('h1tg000041l', '38376304', '43728908', '-'), 3857136, 3878993)]]
    
    # 整体提取序列片段，先选择间隔，再拼装
    not_change_seqs = []
    # 首先检查线段起点到第一个区间的间隙
    not_change_seqs.append([0,merged_final_exchange_regions[0][0][0]]) # [0,2)
    for i in range(len(merged_final_exchange_regions) - 1):
        not_change_seqs.append([merged_final_exchange_regions[i][0][1],merged_final_exchange_regions[i + 1][0][0]]) # [4,5)
    not_change_seqs.append([merged_final_exchange_regions[-1][0][1],targetasm_array_length]) # [7,10) # 统一为真实长度
    
    new_seq = targetasm_fa[not_change_seqs[0][0]:not_change_seqs[0][1]]

    #  record : [updated_sorted_final_match_table,secondasm_fa,secondasm_array_length,secondasm_error_regions]
    
    for i in range(len(merged_final_exchange_regions)):
        record = merged_final_exchange_regions[i][1][0]
        secondasm_fa = contigs_match_table[record][1]
        
        secondasm_block_start = merged_final_exchange_regions[i][1][1]
        secondasm_block_end = merged_final_exchange_regions[i][1][2]
        targetasm_block_start = merged_final_exchange_regions[i][0][0]
        targetasm_block_end = merged_final_exchange_regions[i][0][1]
        new_seq += secondasm_fa[secondasm_block_start:secondasm_block_end]
        new_seq += targetasm_fa[not_change_seqs[i + 1][0]:not_change_seqs[i+1][1]]
    return new_seq,len(new_seq),repair_flag
    

def main():
    
    parser = argparse.ArgumentParser(description="Region repair")
    parser.add_argument("-r", "--repair_file", help="", required=True)
    parser.add_argument("-ts", "--targetasm_fa_file", help="", required=True)
    parser.add_argument("-tf", "--targetasm_error_regions_flag", help="", required=True)
    parser.add_argument("-te", "--targetasm_error_regions_file", help="", required=True)
    
    parser.add_argument("-ssd", "--secondasm_fa_dir", help="", required=True)
    parser.add_argument("-sed", "--secondasm_error_regions_dir", help="", required=True)
    
    parser.add_argument("-omd", "--match_pair_dir",help="", required=True)
    parser.add_argument("-os", "--outfa_file",help="", required=True)
    parser.add_argument("-ol", "--log_file",help="", required=True)
    
    parser.add_argument("-k", "--kmer_size",help="",type=int, default=5000)
    parser.add_argument("-kn", "--kmer_number",help="", type=int,default=2000)
    parser.add_argument("-et", "--extend_length",help="", type=int,default=1000)
    
    args = parser.parse_args()
    
    repair_file = args.repair_file
    targetasm_fa_file = args.targetasm_fa_file
    targetasm_error_regions_flag = args.targetasm_error_regions_flag
    targetasm_error_regions_file = args.targetasm_error_regions_file
    
    secondasm_fa_dir = args.secondasm_fa_dir
    secondasm_error_regions_dir = args.secondasm_error_regions_dir
    match_pair_dir = args.match_pair_dir
    
    
    outfa_file = args.outfa_file
    log_file = args.log_file
    
    kmer_size = args.kmer_size
    kmer_number = args.kmer_number
    extend_length = args.extend_length
    
    
    # repair_file = 'D:/working/AssemblyRepairer/check/region_match.all.xls'
    # targetasm_fa_file = 'D:/working/AssemblyRepairer/check/target_cen_array/chr5_PATERNAL_45671575-58360376.fa'
    # targetasm_error_regions_flag = '+'
    # targetasm_error_regions_file = 'D:/working/AssemblyRepairer/check/target_errors/chr5_PATERNAL_45671575-58360376.xls'
    
    # secondasm_fa_dir = 'D:/working/AssemblyRepairer/check/second_asm_array'
    # secondasm_error_regions_dir = 'D:/working/AssemblyRepairer/check/second_asm_errors'
    
    # match_pair_dir = 'D:/working/AssemblyRepairer/check/out_repair'
    # outfa_file = 'D:/working/AssemblyRepairer/check/out_repair/chr5_PATERNAL_45671575-58360376@+.repaired.fa'
    # log_file = 'D:/working/AssemblyRepairer/check/out_repair/chr5_PATERNAL_45671575-58360376@+.repaired.log'
    
    # kmer_size = 5000
    # kmer_number = 2000
    # extend_length = 1000
    
    if not os.path.exists(match_pair_dir):
        os.mkdir(match_pair_dir)
    
    targetasm_fa = readSeq(targetasm_fa_file,targetasm_error_regions_flag)
    targetasm_array_length = len(targetasm_fa) # # 修改array length为真实array length，不是samtools的end - start，需要额外 + 1
    targetasm_error_regions = []
    error_regions = []
    if targetasm_error_regions_file != '-1':
        targetasm_error_regions = readErrorRegions(targetasm_error_regions_file,targetasm_error_regions_flag,targetasm_array_length)
        error_regions = exchange_regions(targetasm_error_regions,targetasm_array_length,extend_length)
    print(targetasm_error_regions)
    repair_records = []
    
    with open(repair_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            repair_records.append([items[4],items[5],items[6],items[7]])
    
    if len(targetasm_error_regions) == 0:
        return
    
    contigs_match_table = {}
    log_file = open(log_file,'w')
    for i in repair_records:
        # secondasm_fa_file = secondasm_fa_dir + '/' + i[0] + '_' + i[1] + '-' + i[2] + '.fa'
        secondasm_error_regions_flag = i[3]
        secondasm_fa_file = secondasm_fa_dir + '/' + i[0] + ':' + i[1] + '-' + i[2] + '.fa' # linux
        print('---')
        print(targetasm_fa_file)
        print(secondasm_fa_file)
        secondasm_fa = readSeq(secondasm_fa_file,secondasm_error_regions_flag)
        secondasm_array_length = len(secondasm_fa)
        # secondasm_error_regions_file = secondasm_error_regions_dir + '/' + i[0] + '_' + i[1] + '-' + i[2] + '.xls'
        secondasm_error_regions_file = secondasm_error_regions_dir + '/' + i[0] + ':' + i[1] + '-' + i[2] + '.xls'  # linux
        if not os.path.exists(secondasm_error_regions_file):
            secondasm_error_regions_file = '-1'
        secondasm_error_regions = []
        if secondasm_error_regions_file != '-1':
            secondasm_error_regions = readErrorRegions(secondasm_error_regions_file,secondasm_error_regions_flag,secondasm_array_length)
            
        match_file = match_pair_dir + '/' + targetasm_fa_file.split('/')[-1].replace('.fa', '') + '@' + secondasm_fa_file.split('/')[-1].replace('.fa', '')
        if not os.path.exists(match_file):
            sorted_final_match_table = buildKmeMatch(targetasm_fa,secondasm_fa,kmer_size,kmer_number)
            outfile = open(match_file,'w')
            for l in sorted_final_match_table:
                outfile.write(str(l[1][0])+'\t' + str(l[1][1]) + '\t' + l[0] + '\n')
            outfile.close()

        updated_sorted_final_match_table = readKmerPair(match_file)
        
        contigs_match_table[(i[0],i[1],i[2],i[3])] = [updated_sorted_final_match_table,secondasm_fa,secondasm_array_length,secondasm_error_regions]
    print('--------')
    repaired_seq,repaired_seq_len,repair_flag = fixArray_dev(error_regions,contigs_match_table,targetasm_fa,targetasm_array_length,targetasm_error_regions_flag)
    if repair_flag == 0:
        log_file.write('Not repair\n')
    elif repair_flag == 1:
        log_file.write('Repaired All\n')
    elif repair_flag == 3:
        log_file.write('Overlap but repaired Some\n')
    else:
        log_file.write('Boundary limitation\n')
    log_file.close()
    outfa_file = open(outfa_file,'w')
    outfa_file.write('>repaired_' + str(repaired_seq_len) + '_' + targetasm_error_regions_flag+'\n')
    outfa_file.write(repaired_seq+'\n')
    outfa_file.close()
    
    
    
    
    
    
    


if __name__ == '__main__':
    main()
