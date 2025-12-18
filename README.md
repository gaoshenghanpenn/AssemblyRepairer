# AssemblyRepairer
To improve the quality and quantity of the identified centromeres in the HGSVC dataset, we developed an assembly repair method that uses large, unique k-mers (e.g., 5 kbp in length) to identify linked regions between Verkko and Hifiasm contigs and swap the erroneous regions in the Verkko assemblies with the correct, equivalent Hifiasm regions. 

## Dependencies

Python 3.12.0

Development environment: Linux

Development tool: VScode

| Packages           | Version |
| ------------------ | ------- |
| minimap2          | 2.28    |
| pbmm2             | 1.17.0   |
| nucflag              | 0.3.3 |
| samtools         | 1.21   |
| seqkit              | 2.9.0  |

## Quick start
#### Installation

```
Download souce code and install packages using conda
```

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
## Reference
Please cite the following paper when you use AssemblyRepairer in your work

Gao S, Oshima KK, Chuang SC, Loftus M, Montanari A, Gordon DS, Human Genome Structural Variation Consortium, Human Pangenome Reference Consortium, Hsieh P, Konkel MK, Ventura M, Logsdon GA. A global view of human centromere variation and evolution. bioRxiv. 2025. p. 2025.12.09.693231. [doi:10.64898/2025.12.09.693231](https://doi.org/10.64898/2025.12.09.693231)
