# AssemblyRepairer
To improve the quality and quantity of the identified centromeres in the HGSVC dataset, we developed an assembly repair method that uses large, unique k-mers (e.g., 5 kbp in length) to identify linked regions between Verkko and Hifiasm contigs and swap the erroneous regions in the Verkko assemblies with the correct, equivalent Hifiasm regions. 

## Dependencies and Installation
Development environment: Linux  
Development tool: VScode
```
# env.yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python==3.12.0
  - minimap2==2.28
  - pbmm2==1.17.0
  - samtools==1.21
  - seqkit==2.9.0
  - liblzma-devel==5.8.1 # pysam
  - pip
  - pip:
    - nucflag==0.3.3
```
The installation will take approximately a few minutes to ten minutes.  
```
git clone https://github.com/logsdon-lab/AssemblyRepairer.git
conda env create --name AssemblyRepairer -f env.yaml
```

## Quick start
#### Overview
First, run nucflag to detect target assembly errors.

Second, run minimap2 to find the region in second assembly and run nucflag to detect second assembly errors.

Third, based on large, unique k-mers strategy to identify linked regions between two assemblies and swap the erroneous regions.

Fourth, run nucflag to confirm.  

This version is only used for centromere region repairment.
- `target` : Input the coordinates of target regions that needs to be repaired. 

#### Input

* Coordinates file: Regions need to be repaird, like centromere array regions, bed format (contig\tstart\tend\n).
* Target ref file: Assembly that needs repair, like verkko.
* Second assembly file: Second assembly for repair, like hifiasm.
* Hifi reads dir: for nucflag
* Outdir: Output path

```Bash
#For more details, please use python AssemblyRepairer.py target -h
usage: AssemblyRepairer.py target [-h] -c COORDINATES_FILE -taf TARGET_REF_FILE -saf SECOND_ASSEMBLY_FILE -r
                                  HIFI_READS_DIR -o OUTDIR [-rf HIFI_READS_SUFFIX] [-e ERROR_DIR]
                                  [-ig NUCFLAG_IGNORE_REGIONS] [-rb ROLLBACK] [-mq MAPQ_THR]
                                  [-bcm BREAK_CONTIG_MERGE_THR] [-mc MIN_ARRAY_THR] [-om OVERLAP_MERGE_THR]
                                  [-k KMER_SIZE] [-kn KMER_NUMBER] [-et EXTEND_LENGTH] [-t THREADS]

options:
  -h, --help            show this help message and exit
  -c COORDINATES_FILE, --coordinates_file COORDINATES_FILE
                        the coordinates of target regions that needs to be repaired
  -taf TARGET_REF_FILE, --target_ref_file TARGET_REF_FILE
                        Reference assembly needs to be repaired
  -saf SECOND_ASSEMBLY_FILE, --second_assembly_file SECOND_ASSEMBLY_FILE
                        Second assembly use to repair
  -r HIFI_READS_DIR, --hifi_reads_dir HIFI_READS_DIR
                        hifi reads dir
  -o OUTDIR, --outdir OUTDIR
                        outdir
  -rf HIFI_READS_SUFFIX, --hifi_reads_suffix HIFI_READS_SUFFIX
                        hifi reads suffix
  -e ERROR_DIR, --error_dir ERROR_DIR
                        User-defined error regions
  -ig NUCFLAG_IGNORE_REGIONS, --nucflag_ignore_regions NUCFLAG_IGNORE_REGIONS
                        ignore regions for nucflag
  -rb ROLLBACK, --rollback ROLLBACK
                        if region not repair success, then rollback for that error region, 1 use,0 not use
  -mq MAPQ_THR, --MAPQ_thr MAPQ_THR
                        Mapping quality threshold
  -bcm BREAK_CONTIG_MERGE_THR, --break_contig_merge_thr BREAK_CONTIG_MERGE_THR
                        Merge if the interval is smaller than this distance
  -mc MIN_ARRAY_THR, --min_array_thr MIN_ARRAY_THR
                        Minimum array length
  -om OVERLAP_MERGE_THR, --overlap_merge_thr OVERLAP_MERGE_THR
                        merge, If the overlap is less than this value
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        kmer size
  -kn KMER_NUMBER, --kmer_number KMER_NUMBER
                        kmer number for establishing anchor points
  -et EXTEND_LENGTH, --extend_length EXTEND_LENGTH
                        Extend length
  -t THREADS, --threads THREADS
                        Threads for alignment
```
## Demo
Test data is in demo directory, including one HGSVC3 centromere with both verkko and hifiasm assemblies (NA20355 chr16 H2).
```
.
├── coordinate_file.txt 
├── reads # subset hifi reads
│   └── chr16.reads.fa
└── testasm # Considering the running time, this only includes the centromeric region.
    ├── NA20355_chr16_h2tg000041l#1-37317985.fa
    ├── NA20355_chr16_h2tg000041l#1-37317985.fa.fai
    ├── NA20355_chr16_haplotype2-0000183.fa
    └── NA20355_chr16_haplotype2-0000183.fa.fai

```
Running command (The expected demo running time is ~30 minutes, mainly for the pbmm2 alignment)
```
python src/AssemblyRepairer.py target -taf ./demo/testasm/NA20355_chr16_haplotype2-0000183.fa -saf ./demo/testasm/NA20355_chr16_h2tg000041l#1-37317985.fa -c ./demo/coordinate_file.txt -r ./demo/reads -rf fa -o ./out2 -rb 0 -t 24
```
 Output directory, final assembly is repaired_assembly.fa
```
.
├── all.cenarray.fa
├── all.cenarray.fa.paf
├── cenpairs.xls
├── error_cen_array_list.bed
├── NA20355_chr16_haplotype2-0000183.fa
├── NA20355_chr16_haplotype2-0000183.fa.fai
├── nucflag_init_target
│   ├── nucflag.bam
│   ├── nucflag.bam.bai
│   ├── nucflag.misassemblies.bed
│   ├── nucflag.status.bed
│   ├── png
│   │   └── NA20355_chr16_haplotype2-0000183:1-3140550.png
│   └── target_errors
│       └── NA20355_chr16_haplotype2-0000183:1-3140550.xls
├── nucflag_repaired_asm
│   ├── nucflag.bam
│   ├── nucflag.bam.bai
│   ├── nucflag.misassemblies.bed
│   ├── nucflag.status.bed
│   ├── png
│   │   └── NA20355_chr16_haplotype2-0000183:1-3145645.png
│   └── repaired_asm_errors
├── nucflag_second_asm
│   ├── nucflag.bam
│   ├── nucflag.bam.bai
│   ├── nucflag.misassemblies.bed
│   ├── nucflag.status.bed
│   ├── png
│   │   └── NA20355_chr16_h2tg000041l#1-37317985:1-2861716.png
│   └── second_asm_errors
│       └── NA20355_chr16_h2tg000041l#1-37317985:1-2861716.xls
├── out_repair
│   ├── NA20355_chr16_haplotype2-0000183:1-3140550@NA20355_chr16_h2tg000041l#1-37317985:1-2861716
│   ├── NA20355_chr16_haplotype2-0000183:1-3140550@+.repaired.fa
│   ├── NA20355_chr16_haplotype2-0000183:1-3140550@+.repaired.log
│   └── NA20355_chr16_haplotype2-0000183:1-3140550@+.repair.xls
├── region_match.all.xls
├── repaired_assembly.fa # final repaired fa file
├── repaired.bed # region update between after and before.
├── seconda_asm_array.nucflag.bed
├── second_asm_array
│   └── NA20355_chr16_h2tg000041l#1-37317985:1-2861716.fa
└── target_cen_array
    └── NA20355_chr16_haplotype2-0000183:1-3140550.fa
```
Before repair
![before](img/NA20355_chr16_haplotype2-0000183_1-3140550.png)  
hifiasm array
![hifiasm](img/NA20355_chr16_h2tg000041l#1-37317985_1-2861716.png)  
After repair
![hifiasm](img/NA20355_chr16_haplotype2-0000183_1-3145645.png)  

## Citing 
**Gao S, Oshima KK**, Chuang SC, Loftus M, Montanari A, Gordon DS, Human Genome Structural Variation Consortium, Human Pangenome Reference Consortium, Hsieh P, Konkel MK, Ventura M, Logsdon GA. A global view of human centromere variation and evolution. bioRxiv. 2025. p. 2025.12.09.693231. [doi:10.64898/2025.12.09.693231](https://doi.org/10.64898/2025.12.09.693231)