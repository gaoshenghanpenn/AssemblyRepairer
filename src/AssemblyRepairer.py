import argparse
import os
from repair_utils import *

script_path = os.path.split(os.path.realpath(__file__))[0]

def targetRegionRepair(args):
    coordinates_file = args.coordinates_file
    target_ref_file = args.target_ref_file
    second_assembly_file = args.second_assembly_file
    hifi_reads_dir = args.hifi_reads_dir
    
    outdir = args.outdir
    
    rollback = args.rollback
    
    hifi_reads_suffix = args.hifi_reads_suffix
    
    error_dir = args.error_dir
    nucflag_ignore_regions = args.nucflag_ignore_regions
    
    MAPQ_thr=args.MAPQ_thr
    break_contig_merge_thr=args.break_contig_merge_thr
    min_array_thr=args.min_array_thr
    overlap_merge_thr=args.overlap_merge_thr
    
    kmer_size=args.kmer_size
    kmer_number=args.kmer_number
    extend_length=args.extend_length
    
    mapping_threads = args.threads
    hifi_reads_mapping_threads = args.threads
    
    
    
    # coordinates_file = '/project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/AssemblyRepairer/testData/cen_array_list.bed'
    # target_ref_file = '/project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/AssemblyRepairer/testData/HG00096-asm-renamed-reort.fa'
    # second_assembly_file = '/project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/AssemblyRepairer/testData/HG00096-asm-renamed-reort.hifiasm.fa'
    # hifi_reads_dir = '/project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/data/hifi/HG00096'
    # hifi_reads_suffix = '.fastq.gz'
    
    # error_dir = ''
    # nucflag_ignore_regions = '/project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/AssemblyRepairer/testData/HG00096_correct_ALR_regions.rm.simple.bed'
    
    # outdir = '/project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/AssemblyRepairer/testData/outdir'
    
    
    # MAPQ_thr=20
    # break_contig_merge_thr=5000
    # min_array_thr=100000
    # overlap_merge_thr=0.3
    
    # kmer_size=5000
    # kmer_number=2000
    # extend_length=1000
    
    # mapping_threads = 10
    # hifi_reads_mapping_threads = 25
    
    nucflag_config = script_path + '/nucflag.toml'
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    all_cen_arrays = []
    with open(coordinates_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            array_name = items[0] + ':' + items[1] + '-' + items[2]
            all_cen_arrays.append(array_name)
    
    # 0. run nucflag for coordinates_file, get error regions
    # or input error bed
    if error_dir == '':
        # run nucflag
        files = os.listdir(hifi_reads_dir)
        hifi_reads_files = []
        for i in files:
            if i.endswith(hifi_reads_suffix):
                hifi_reads_files.append(hifi_reads_dir + '/' + i)
        out_nucflag_dir = outdir + '/nucflag_init_target'
        if not os.path.exists(out_nucflag_dir):
            os.mkdir(out_nucflag_dir)
        error_dir = out_nucflag_dir + '/target_errors'
        run_nucflag(hifi_reads_files,target_ref_file,out_nucflag_dir,hifi_reads_mapping_threads,nucflag_ignore_regions,nucflag_config,coordinates_file,error_dir)
    
    # build error array file
    error_cen_arrays = []
    error_coordinates_file = outdir + '/error_cen_array_list.bed'
    files = os.listdir(error_dir)
    out_error_coordinates_file = open(error_coordinates_file,'w')
    for i in files:
        array = i.replace('.xls', '')
        error_cen_arrays.append(array)
        info = array.split(':')
        contig = info[0]
        region = info[1].split('-')
        out_error_coordinates_file.write(contig + '\t' + region[0] + '\t' + region[1] + '\n')
    out_error_coordinates_file.close()
        
    
    # 1. get target cen array
    target_cen_array = outdir + '/target_cen_array'
    if not os.path.exists(target_cen_array):
        os.mkdir(target_cen_array)
    run_getTargetArray(target_ref_file,error_coordinates_file,outdir,target_cen_array)
    
    # 2. mapping cen array to second assembly
    cen_array_fa_file = outdir + '/all.cenarray.fa'
    paf_file = cen_array_fa_file + '.paf'
    cmd = 'minimap2 -x asm5 -t ' + str(mapping_threads) +' --eqx --cs --secondary=no -s 25000 -K 8G ' + second_assembly_file + ' ' + cen_array_fa_file + ' > ' + paf_file
    os.system(cmd)
    
    # 3. get matching records from paf
    matching_records_file = outdir + '/cenpairs.xls'
    getFilterRecord(paf_file,error_coordinates_file,matching_records_file,MAPQ_thr,break_contig_merge_thr,min_array_thr,overlap_merge_thr)
    
    seconda_asm_array_bed = outdir + '/seconda_asm_array.nucflag.bed'
    repair_table = outdir + '/region_match.all.xls'
    buildInitNucflagRegions(matching_records_file,seconda_asm_array_bed,repair_table)
    
    # build second asm array fasta
    second_asm_array = outdir + '/second_asm_array' 
    if not os.path.exists(second_asm_array):
        os.mkdir(second_asm_array)
    with open(seconda_asm_array_bed,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            cmd = 'samtools faidx ' + second_assembly_file + ' ' + items[0] + ':' + items[1] + '-' + items[2] + ' > ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa'
            os.system(cmd)
            cmd = 'seqkit seq -w 0 ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa' + ' > ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.fa'
            os.system(cmd)
            cmd = 'rm ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa'
            os.system(cmd)
    
    # 4. nucflag check for second assembly 
    out_nucflag_dir = outdir + '/nucflag_second_asm'
    if not os.path.exists(out_nucflag_dir):
        os.mkdir(out_nucflag_dir)
    
    second_error_dir = out_nucflag_dir + '/second_asm_errors'
    run_nucflag(hifi_reads_files,second_assembly_file,out_nucflag_dir,hifi_reads_mapping_threads,'',nucflag_config,seconda_asm_array_bed,second_error_dir)
    
    # 5. run array repair
    repairdir = outdir + '/out_repair'
    if not os.path.exists(repairdir):
        os.mkdir(repairdir)
    print('start repair')
    script = script_path + '/region_repairment.py' 
    run_repair(script, repair_table,target_cen_array,second_asm_array,error_dir,second_error_dir,repairdir,kmer_size,kmer_number,extend_length)
    
    # 6. run assembly repair
    script = script_path + '/repair_assembly.py' 
    outbed = outdir + '/repaired.bed'
    repaired_assembly = outdir + '/repaired_assembly.fa'
    cmd = 'python ' + script + ' -r ' + target_ref_file + ' -a ' + repairdir + ' -o ' + outdir + ' -or ' + outbed + ' -of ' + repaired_assembly
    os.system(cmd)
    
    
    # 7. run final nucflag check
    out_nucflag_dir = outdir + '/nucflag_repaired_asm'
    if not os.path.exists(out_nucflag_dir):
        os.mkdir(out_nucflag_dir)
    final_error_dir = out_nucflag_dir + '/repaired_asm_errors'
    run_nucflag(hifi_reads_files,repaired_assembly,out_nucflag_dir,hifi_reads_mapping_threads,'',nucflag_config,outbed,final_error_dir)
    
    # 8. rollback
    if rollback == 1:
        repaired_regions = {}
        repaired_state = {}
        with open(outbed,'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')
                new_region = items[0] + ':' + items[1] + '-' + items[2]
                old_region = items[0] + ':' + items[4] + '-' + items[5]
                repaired_regions[new_region] = old_region
        
        error_arrays = set()
        with open(out_nucflag_dir + '/nucflag.misassemblies.bed', 'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                if 'HET' in line:
                    continue
                items = line.split('\t')[0]
                error_arrays.add(items)
        
        
        for i in repaired_regions.keys():
            if i in error_arrays:
                repaired_state[i] = ['error',repaired_regions[i]]
            else:
                repaired_state[i] = ['success',repaired_regions[i]]
        
        outRollback_file_name = outdir + '/rollback.bed'
        outRollback_file = open(outRollback_file_name,'w')
        for i in repaired_state.keys():
            outRollback_file.write(i + '\t' + repaired_state[i][1] + '\t' + repaired_state[i][0] + '\n' )
        
        outRollback_file.close()        
                
        # cmd rollback
        script = script_path + '/rollback.py' 
        repaired_assembly = outdir + '/repaired_assembly.rollback.fa'
        outbed = outdir + '/repaired.rollback.bed'
        cmd = 'python ' + script + ' -rb ' + outRollback_file_name + ' -r ' + target_ref_file + ' -a ' + repairdir + ' -o ' + outdir + ' -or ' + outbed + ' -of ' + repaired_assembly
        os.system(cmd)


def repairOnly(args):
    
    repair_table = args.repair_table
    target_cen_array = args.target_cen_array
    second_asm_array = args.second_asm_array
    error_dir = args.error_dir
    second_error_dir = args.second_error_dir
    target_ref_file = args.target_ref_file
    outdir = args.outdir
    
    kmer_size=args.kmer_size
    kmer_number=args.kmer_number
    extend_length=args.extend_length
    
    
    
    # 5. run array repair
    repairdir = outdir + '/out_repair'
    if not os.path.exists(repairdir):
        os.mkdir(repairdir)
    print('start repair')
    script = script_path + '/region_repairment.py' 
    run_repair(script, repair_table,target_cen_array,second_asm_array,error_dir,second_error_dir,repairdir,kmer_size,kmer_number,extend_length)
    
    # 6. run assembly repair
    script = script_path + '/repair_assembly.py' 
    outbed = outdir + '/repaired.bed'
    repaired_assembly = outdir + '/repaired_assembly.fa'
    cmd = 'python ' + script + ' -r ' + target_ref_file + ' -a ' + repairdir + ' -o ' + outdir + ' -or ' + outbed + ' -of ' + repaired_assembly
    os.system(cmd)
    
def matchBuilder(args):
    
    cen_array_fa_file = args.cen_array_fa_file
    second_assembly_file = args.second_assembly_file
    error_region_list_file = args.error_region_list_file
    
    outdir = args.outdir
    
    MAPQ_thr=args.MAPQ_thr
    break_contig_merge_thr=args.break_contig_merge_thr
    min_array_thr=args.min_array_thr
    overlap_merge_thr=args.overlap_merge_thr
    
    
    mapping_threads = args.threads
    
    # 2. mapping cen array to second assembly
    paf_file = cen_array_fa_file + '.paf'
    cmd = 'minimap2 -x asm5 -t ' + str(mapping_threads) +' --eqx --cs --secondary=no -s 25000 -K 8G ' + second_assembly_file + ' ' + cen_array_fa_file + ' > ' + paf_file
    os.system(cmd)
    
    # 3. get matching records from paf
    matching_records_file = outdir + '/cenpairs.xls'
    getFilterRecord(paf_file,error_region_list_file,matching_records_file,MAPQ_thr,break_contig_merge_thr,min_array_thr,overlap_merge_thr)
    
    seconda_asm_array_bed = outdir + '/seconda_asm_array.nucflag.bed'
    repair_table = outdir + '/region_match.all.xls'
    buildInitNucflagRegions(matching_records_file,seconda_asm_array_bed,repair_table)


def rollBack(args):
    
    rollback_file = args.rollback_file
    target_ref_file = args.target_ref_file
    repairdir = args.repairdir
    
    outdir = args.outdir
    # cmd rollback
    script = script_path + '/rollback.py' 
    repaired_assembly = outdir + '/repaired_assembly.rollback.fa'
    outbed = outdir + '/repaired.rollback.bed'
    cmd = 'python ' + script + ' -rb ' + rollback_file + ' -r ' + target_ref_file + ' -a ' + repairdir + ' -o ' + outdir + ' -or ' + outbed + ' -of ' + repaired_assembly
    os.system(cmd)


def whole_genome(args):
    # end to end pipline: only support nucflag check
    
    target_ref_file = args.target_ref_file
    second_assembly_file = args.second_assembly_file
    hifi_reads_dir = args.hifi_reads_dir
    
    outdir = args.outdir
    
    rollback = args.rollback
    
    hifi_reads_suffix = args.hifi_reads_suffix
    
    error_dir = args.error_dir
    
    nucflag_ignore_regions = args.nucflag_ignore_regions
    
    merge_distance = args.merge_distance
    
    MAPQ_thr=args.MAPQ_thr
    break_contig_merge_thr=args.break_contig_merge_thr
    min_array_thr=args.min_array_thr
    overlap_merge_thr=args.overlap_merge_thr
    
    kmer_size=args.kmer_size
    kmer_number=args.kmer_number
    extend_length=args.extend_length
    
    mapping_threads = args.threads
    hifi_reads_mapping_threads = args.threads
    
    
    
    nucflag_config = script_path + '/nucflag.toml'
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    files = os.listdir(hifi_reads_dir)
    hifi_reads_files = []
    for i in files:
        if i.endswith(hifi_reads_suffix):
            hifi_reads_files.append(hifi_reads_dir + '/' + i)
    out_nucflag_dir = outdir + '/nucflag_init_target'
    if not os.path.exists(out_nucflag_dir):
        os.mkdir(out_nucflag_dir)
    error_dir = out_nucflag_dir + '/target_errors_0'
    # run_nucflag(hifi_reads_files,target_ref_file,out_nucflag_dir,hifi_reads_mapping_threads,nucflag_ignore_regions,nucflag_config,'',error_dir)
    
    # 构建错误区间，然后整理生成target_errors
    error_file = out_nucflag_dir + '/nucflag.misassemblies.bed'
    cmd = 'samtools faidx ' + target_ref_file
    os.system(cmd)
    fai_file = target_ref_file + '.fai'
    error_region_file = outdir + '/error_regions.bed'
    buildAllErrorRegions(error_file,fai_file,error_region_file,merge_distance)
    
    error_regions = {}
    with open(error_region_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            error_regions[(items[0],int(items[1]), int(items[2]))] = []
    
    all_error_regions = []
    with open(error_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            if 'HET' in line:
                continue
            items = line.split('\t')
            all_error_regions.append([items[0],int(items[1]) + 1,int(items[2])])
    
    for i in all_error_regions:
        for j in error_regions.keys():
            if i[0] == j[0]:
                if i[1] >= j[1] and i[1] <= j[2]:
                    error_regions[j].append([i[1],i[2]])
                    break
    
    error_dir = out_nucflag_dir + '/target_errors'
    if not os.path.exists(error_dir):
        os.mkdir(error_dir)
        
    for i in error_regions.keys():
        error_file = error_dir + '/' + i[0] + ':' + str(i[1]) + '-' + str(i[2]) + '.xls'
        error_file = open(error_file,'w')
        for j in error_regions[i]:
            start = j[0] - i[1]
            end = j[1] - i[1]
            error_file.write(str(start) + '\t'  +str(end) + '\n')
        error_file.close()
    
    
    # build error array file
    error_cen_arrays = []
    error_coordinates_file = outdir + '/error_cen_array_list.bed'
    files = os.listdir(error_dir)
    out_error_coordinates_file = open(error_coordinates_file,'w')
    for i in files:
        array = i.replace('.xls', '')
        error_cen_arrays.append(array)
        info = array.split(':')
        contig = info[0]
        region = info[1].split('-')
        out_error_coordinates_file.write(contig + '\t' + region[0] + '\t' + region[1] + '\n')
    out_error_coordinates_file.close()
        
    
    # 1. get target cen array
    target_cen_array = outdir + '/target_cen_array'
    if not os.path.exists(target_cen_array):
        os.mkdir(target_cen_array)
    run_getTargetArray(target_ref_file,error_coordinates_file,outdir,target_cen_array)
    
    # 2. mapping cen array to second assembly
    cen_array_fa_file = outdir + '/all.cenarray.fa'
    paf_file = cen_array_fa_file + '.paf'
    cmd = 'minimap2 -x asm5 -t ' + str(mapping_threads) +' --eqx --cs --secondary=no -s 25000 -K 8G ' + second_assembly_file + ' ' + cen_array_fa_file + ' > ' + paf_file
    os.system(cmd)
    
    # 3. get matching records from paf
    matching_records_file = outdir + '/cenpairs.xls'
    getFilterRecord(paf_file,error_coordinates_file,matching_records_file,MAPQ_thr,break_contig_merge_thr,min_array_thr,overlap_merge_thr)
    
    seconda_asm_array_bed = outdir + '/seconda_asm_array.nucflag.bed'
    repair_table = outdir + '/region_match.all.xls'
    buildInitNucflagRegions(matching_records_file,seconda_asm_array_bed,repair_table)
    
    # build second asm array fasta
    second_asm_array = outdir + '/second_asm_array' 
    if not os.path.exists(second_asm_array):
        os.mkdir(second_asm_array)
    with open(seconda_asm_array_bed,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            cmd = 'samtools faidx ' + second_assembly_file + ' ' + items[0] + ':' + items[1] + '-' + items[2] + ' > ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa'
            os.system(cmd)
            cmd = 'seqkit seq -w 0 ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa' + ' > ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.fa'
            os.system(cmd)
            cmd = 'rm ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa'
            os.system(cmd)
    
    # 4. nucflag check for second assembly 
    out_nucflag_dir = outdir + '/nucflag_second_asm'
    if not os.path.exists(out_nucflag_dir):
        os.mkdir(out_nucflag_dir)
    
    second_error_dir = out_nucflag_dir + '/second_asm_errors'
    run_nucflag(hifi_reads_files,second_assembly_file,out_nucflag_dir,hifi_reads_mapping_threads,'',nucflag_config,seconda_asm_array_bed,second_error_dir)
    
    # 5. run array repair
    repairdir = outdir + '/out_repair'
    if not os.path.exists(repairdir):
        os.mkdir(repairdir)
    print('start repair')
    script = script_path + '/region_repairment.py' 
    run_repair(script, repair_table,target_cen_array,second_asm_array,error_dir,second_error_dir,repairdir,kmer_size,kmer_number,extend_length)
    
    # 6. run assembly repair
    script = script_path + '/repair_assembly.py' 
    outbed = outdir + '/repaired.bed'
    repaired_assembly = outdir + '/repaired_assembly.fa'
    cmd = 'python ' + script + ' -r ' + target_ref_file + ' -a ' + repairdir + ' -o ' + outdir + ' -or ' + outbed + ' -of ' + repaired_assembly
    os.system(cmd)
    
    
    # 7. run final nucflag check
    out_nucflag_dir = outdir + '/nucflag_repaired_asm'
    if not os.path.exists(out_nucflag_dir):
        os.mkdir(out_nucflag_dir)
    final_error_dir = out_nucflag_dir + '/repaired_asm_errors'
    run_nucflag(hifi_reads_files,repaired_assembly,out_nucflag_dir,hifi_reads_mapping_threads,'',nucflag_config,outbed,final_error_dir)
    
    # 8. rollback
    if rollback == 1:
        repaired_regions = {}
        repaired_state = {}
        with open(outbed,'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')
                new_region = items[0] + ':' + items[1] + '-' + items[2]
                old_region = items[0] + ':' + items[4] + '-' + items[5]
                repaired_regions[new_region] = old_region
        
        error_arrays = set()
        with open(out_nucflag_dir + '/nucflag.misassemblies.bed', 'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                if 'HET' in line:
                    continue
                items = line.split('\t')[0]
                error_arrays.add(items)
        
        
        for i in repaired_regions.keys():
            if i in error_arrays:
                repaired_state[i] = ['error',repaired_regions[i]]
            else:
                repaired_state[i] = ['success',repaired_regions[i]]
        
        outRollback_file_name = outdir + '/rollback.bed'
        outRollback_file = open(outRollback_file_name,'w')
        for i in repaired_state.keys():
            outRollback_file.write(i + '\t' + repaired_state[i][1] + '\t' + repaired_state[i][0] + '\n' )
        
        outRollback_file.close()        
                
        # cmd rollback
        script = script_path + '/rollback.py' 
        repaired_assembly = outdir + '/repaired_assembly.final.fa'
        outbed = outdir + '/repaired.final.bed'
        cmd = 'python ' + script + ' -rb ' + outRollback_file_name + ' -r ' + target_ref_file + ' -a ' + repairdir + ' -o ' + outdir + ' -or ' + outbed + ' -of ' + repaired_assembly
        os.system(cmd)
    
    
    
    
    
    


def main():
    parser = argparse.ArgumentParser(description="AssemblyRepairer: Automated assembly repair tool")

    subparsers = parser.add_subparsers(help="Modes: target, whole_genome (Under development), match_builder,repair_only, rollback")
    
    # python AssemblyRapairer.py target
    parser_target = subparsers.add_parser('target', help='repair target regions in assembly')
    parser_target.add_argument("-c", "--coordinates_file", help="the coordinates of target regions that needs to be repaired", required=True)
    parser_target.add_argument("-taf", "--target_ref_file", help="Reference assembly needs to be repaired", required=True)
    parser_target.add_argument("-saf", "--second_assembly_file", help="Second assembly use to repair", required=True)
    parser_target.add_argument("-r", "--hifi_reads_dir", help="hifi reads dir", required=True)
    parser_target.add_argument("-o", "--outdir", help="outdir", required=True)
    parser_target.add_argument("-rf", "--hifi_reads_suffix", help="hifi reads suffix", default='.fastq.gz')
    parser_target.add_argument("-e", "--error_dir", help="User-defined error regions", default='')
    parser_target.add_argument("-ig", "--nucflag_ignore_regions", help="ignore regions for nucflag", default='')
    parser_target.add_argument("-rb", "--rollback", help="if region not repair success, then rollback for that error region, 1 use,0 not use",type=int, default=1)
    parser_target.add_argument("-mq", "--MAPQ_thr",help="Mapping quality threshold",type=int, default=20)
    parser_target.add_argument("-bcm", "--break_contig_merge_thr",help="Merge if the interval is smaller than this distance", type=int,default=5000)
    parser_target.add_argument("-mc", "--min_array_thr",help="Minimum array length", type=int,default=100000)
    parser_target.add_argument("-om", "--overlap_merge_thr",help="merge, If the overlap is less than this value", type=float,default=0.3)
    
    parser_target.add_argument("-k", "--kmer_size",help="kmer size",type=int, default=5000)
    parser_target.add_argument("-kn", "--kmer_number",help="kmer number for establishing anchor points", type=int,default=2000)
    parser_target.add_argument("-et", "--extend_length",help="Extend length", type=int,default=1000)
    parser_target.add_argument("-t", "--threads",help="Threads for alignment", type=int,default=1)
    parser_target.set_defaults(func=targetRegionRepair)
    
    
    # python AssemblyRapairer.py repair_only    
    parser_repair_only = subparsers.add_parser('repair_only', help='Minimum Operation to repair target regions in assembly')
    
    parser_repair_only.add_argument("-r", "--repair_table", help="the records of repairment", required=True)
    parser_repair_only.add_argument("-taf", "--target_ref_file", help="Reference assembly needs to be repaired", required=True)
    
    parser_repair_only.add_argument("-tca", "--target_cen_array", help="the directory contain region fasta file of target assemblies", required=True)
    parser_repair_only.add_argument("-saa", "--second_asm_array", help="the directory contain region fasta file of second assemblies", required=True)
    
    parser_repair_only.add_argument("-er", "--error_dir", help="the directory contain relative coordinates of the error in target regions", required=True)
    parser_repair_only.add_argument("-ser", "--second_error_dir", help="the directory contain relative coordinates of the error in second regions", required=True)
    
    parser_repair_only.add_argument("-o", "--outdir", help="outdir", required=True)
    
    parser_repair_only.add_argument("-k", "--kmer_size",help="kmer size",type=int, default=5000)
    parser_repair_only.add_argument("-kn", "--kmer_number",help="kmer number for establishing anchor points", type=int,default=2000)
    parser_repair_only.add_argument("-et", "--extend_length",help="Extend length", type=int,default=1000)
    
    parser_repair_only.set_defaults(func=repairOnly)
    
    # python AssemblyRapairer.py match_builder  
    parser_match_builder = subparsers.add_parser('match_builder', help='Used to build match table between two assemblies')
    parser_match_builder.add_argument("-cf", "--cen_array_fa_file", help="A fasta file contains all sequences", required=True)
    parser_match_builder.add_argument("-saf", "--second_assembly_file", help="Second assembly use to repair", required=True)
    parser_match_builder.add_argument("-erl", "--error_region_list_file", help="A file contains all target regions, like chr:start-end", required=True)
    
    parser_match_builder.add_argument("-o", "--outdir", help="outdir", required=True)
    
    parser_match_builder.add_argument("-mq", "--MAPQ_thr",help="Mapping quality threshold",type=int, default=20)
    parser_match_builder.add_argument("-bcm", "--break_contig_merge_thr",help="Merge if the interval is smaller than this distance", type=int,default=5000)
    parser_match_builder.add_argument("-mc", "--min_array_thr",help="Minimum array length", type=int,default=100000)
    parser_match_builder.add_argument("-om", "--overlap_merge_thr",help="merge, If the overlap is less than this value", type=float,default=0.3)
    
    parser_match_builder.add_argument("-t", "--mapping_threads",help="Threads for alignment", type=int,default=1)
    parser_match_builder.set_defaults(func=matchBuilder)
    
    
    # python AssemblyRapairer.py rollback  
    parser_rollback = subparsers.add_parser('rollback', help='Used to rollback error repairment')
    parser_rollback.add_argument("-rb", "--rollback_file", help="A config file for rollback", required=True)
    parser_rollback.add_argument("-taf", "--target_ref_file", help="Reference assembly needs to be repaired", required=True)
    parser_rollback.add_argument("-ro", "--repairdir", help="repair output directory, which contains repaired arrays", required=True)
    parser_rollback.add_argument("-o", "--outdir", help="outdir", required=True)
    parser_rollback.set_defaults(func=rollBack)
    
    # python AssemblyRapairer.py whole_genome     
    parser_whole_genome = subparsers.add_parser('whole_genome', help='repair assembly')
    parser_whole_genome.add_argument("-taf", "--target_ref_file", help="Reference assembly needs to be repaired", required=True)
    parser_whole_genome.add_argument("-saf", "--second_assembly_file", help="Second assembly use to repair", required=True)
    parser_whole_genome.add_argument("-r", "--hifi_reads_dir", help="hifi reads dir", required=True)
    parser_whole_genome.add_argument("-o", "--outdir", help="outdir", required=True)
    parser_whole_genome.add_argument("-rf", "--hifi_reads_suffix", help="hifi reads suffix", default='.fastq.gz')
    parser_whole_genome.add_argument("-e", "--error_dir", help="User-defined error regions", default='')
    parser_whole_genome.add_argument("-ig", "--nucflag_ignore_regions", help="ignore regions for nucflag", default='')
    
    parser_whole_genome.add_argument("-md", "--merge_distance", help="error regions merge distance and extend half, default 500k",type=int, default=500000)
    
    parser_whole_genome.add_argument("-rb", "--rollback", help="if region not repair success, then rollback for that error region, 1 use,0 not use",type=int, default=1)
    parser_whole_genome.add_argument("-mq", "--MAPQ_thr",help="Mapping quality threshold",type=int, default=20)
    parser_whole_genome.add_argument("-bcm", "--break_contig_merge_thr",help="Merge if the interval is smaller than this distance", type=int,default=5000)
    parser_whole_genome.add_argument("-mc", "--min_array_thr",help="Minimum array length", type=int,default=100000)
    parser_whole_genome.add_argument("-om", "--overlap_merge_thr",help="merge, If the overlap is less than this value", type=float,default=0.3)
    
    parser_whole_genome.add_argument("-k", "--kmer_size",help="kmer size",type=int, default=5000)
    parser_whole_genome.add_argument("-kn", "--kmer_number",help="kmer number for establishing anchor points", type=int,default=2000)
    parser_whole_genome.add_argument("-et", "--extend_length",help="Extend length", type=int,default=1000)
    parser_whole_genome.add_argument("-t", "--threads",help="Threads for alignment", type=int,default=1)
    parser_whole_genome.set_defaults(func=whole_genome)
    
    
    args = parser.parse_args()
    args.func(args)
    
    
    
    

if __name__ == '__main__':
    main()