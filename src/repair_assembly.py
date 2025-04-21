
import os
import argparse

def get_reverse_complement(dna_sequence):
    # 定义碱基互补对
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
    # 生成互补链
    complement_sequence = ''.join(complement[base] for base in dna_sequence)
    # 反转互补链得到反链
    reverse_complement_sequence = complement_sequence[::-1]
    return reverse_complement_sequence

def read_ori_assembly(file):
    fa = {}
    header = ''
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            header = line[1:]
            seq = f.readline()[:-1]
            fa[header] = seq
    return fa


def main():
    
    parser = argparse.ArgumentParser(description="Assembly repair")
    parser.add_argument("-r", "--ori_assembly_file", help="", required=True)
    parser.add_argument("-a", "--array_dir", help="", required=True)
    parser.add_argument("-o", "--outdir", help="", required=True)
    parser.add_argument("-or", "--out_array_region_file", help="", required=True)
    parser.add_argument("-of", "--out_file", help="", required=True)
    
    args = parser.parse_args()
    
    ori_assembly_file = args.ori_assembly_file
    array_dir = args.array_dir
    outdir = args.outdir
    out_array_region_file = args.out_array_region_file
    outfile = args.out_file
    
    # build 1 line assembly
    assembly_name = ori_assembly_file.split('/')[-1]
    cmd = 'seqkit seq -w 0 ' + ori_assembly_file + '  >  ' + outdir + '/' + assembly_name
    os.system(cmd)
    cmd = 'samtools faidx ' +  outdir + '/' + assembly_name
    os.system(cmd)
    
    # assembly_name = 'HG00171.fa'
    # outdir = 'D:/working/AssemblyRepairer/test_assembly_repair'
    # array_dir = 'D:/working/AssemblyRepairer/test_assembly_repair/out_repair'
    # ori_assembly_file =  outdir + '/' + assembly_name
    # out_array_region_file = outdir + '/repaired.bed' 
    # outfile = outdir + '/repaired.fa' 
    
    fa_seqs = read_ori_assembly(outdir + '/' + assembly_name)
    print(ori_assembly_file)
    files = os.listdir(array_dir)
    arrays = {}
    for i in files:
        if i.endswith('.fa'):
            array_name = i.split('.')[0].split('@')[0]
            contig_name = array_name.split(':')[0] # linux
            # contig_name = array_name.split('_')[0] + '_' + array_name.split('_')[1] + '_' + array_name.split('_')[2]  
            if contig_name not in arrays.keys():
                arrays[contig_name] = []
            print(array_name)
            region = array_name.split(':')[1].split('-') # linux str(), [start,end]
            # region = array_name.split('_')[-1].split('-') 
            region = [int(region[0]),int(region[1])]
            fa_file = array_dir + '/' + i
            with open(fa_file,'r') as f:
                line = f.readline()[:-1]
                header = line[1:]
                flag = header.split('_')[-1]
                seq = f.readline()[:-1]
                # 由于一个contig可能存在多个区间，修改为list保存
                # 合并arrays和arrays_info，便于排序
                arrays[contig_name].append([region,flag,seq])
    
    sorted_arrays = {}
    for i in arrays.keys():
        sorted_arrays[i] = sorted(arrays[i],key=lambda x:x[0][0])
    
    arrays = sorted_arrays
    sum_base_number = 0
    outfile = open(outfile,'w')
    out_array_region_file = open(out_array_region_file,'w')
    for i in fa_seqs.keys():
        if i not in arrays.keys():
            outfile.write('>' + i+'\n')
            outfile.write(fa_seqs[i] + '\n')
            sum_base_number += len(fa_seqs[i])
        else:
            # 变换chr
            arrays_in_one_contig = arrays[i]
            # 首先遍历获得拼接的ref区间
            start = 0
            regions_in_ref = []
            for j in arrays_in_one_contig:
                region = j[0]
                end = region[0] - 1
                regions_in_ref.append([start,end])
                start = region[1]
            regions_in_ref.append([start,len(fa_seqs[i])])
            # 存在，在最开始就是SD的区间[0,0]
            # 遍历拼接，先拼接一个ref，在拼接一个target + ref
            new_regions = []
            init_seq_region = regions_in_ref[0]
            new_seq = fa_seqs[i][init_seq_region[0] : init_seq_region[1]]
            for j in range(len(arrays_in_one_contig)):
                new_start = len(new_seq) + 1
                new_end = len(new_seq) + len(arrays_in_one_contig[j][2]) 
                new_seq += arrays_in_one_contig[j][2]

                ref_region = regions_in_ref[j + 1]
                new_seq += fa_seqs[i][ref_region[0] : ref_region[1]]
                region = arrays_in_one_contig[j][0]
                new_regions.append([new_start,new_end,new_end - new_start,region[0],region[1]]) # samtools region

            outfile.write('>' + i + '\n')
            outfile.write(new_seq + '\n')
            sum_base_number += len(new_seq)
            for j in new_regions:
                out_array_region_file.write(i+'\t' + str(j[0]) + '\t' + str(j[1]) + '\t' + str(j[2]) + '\t' + str(j[3]) + '\t' + str(j[4]) + '\n')
    
    outfile.close()
    out_array_region_file.close()
    


if __name__ == '__main__':
    main()
