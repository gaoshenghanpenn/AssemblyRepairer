# AssemblyRepairer
To improve the quality and quantity of the identified centromeres in the HGSVC dataset, we developed an assembly repair method that uses large, unique k-mers (e.g., 5 kbp in length) to identify linked regions between Verkko and Hifiasm contigs and swap the erroneous regions in the Verkko assemblies with the correct, equivalent Hifiasm regions. 

## Dependencies

Python 3.12.0

Development environment: Linux

Development tool: VScode

| Packages           | Version |
| ------------------ | ------- |
| minimap2          | 2.28    |
| pbmm2             | 1.14.99   |
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

For detail usage, run `python AssemblyRepairer.py -h`

## Contact

