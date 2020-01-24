# GATK* pipeline for Allelic Counting
([GATK](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)-like) data preprocessing for Allele-specific expression analysis on RNA-seq data.

This pipeline was used as the main procedure for creating tables of gene allelic counts for the further AI analysis with [QCumber](https://github.com/gimelbrantlab/QCumber) in paper

>_"Unexpected variability of allelic imbalance estimates from RNA sequencing", Mendelevich A.*, Vinogradova S.*, Gupta S., Mironov A., Synyaev S., Gimelbrant A._

The pipeline can be divided into two prime parts:

1. **Reference preparation**
  constructs individual paternal and maternal genomes, creates heterozigous VCF.

2. **Creation of tables with Allelic Counts**
  maps the reads from RNA-seq experiments to these genomes and counts the number of reads which map to either the reference or alternate allele at each heterozygous SNP, estimates allelic imbalance for individual genes summarizing information from SNPs.

Please find **manuals / worked examples** at **[Wiki page](https://github.com/gimelbrantlab/GATKstar/wiki)** of this repository.

The resulting allelic counts tables can be used for allelic imbalance analysis via [QCumber](https://github.com/gimelbrantlab/QCumber).

![pic](https://github.com/gimelbrantlab/GATKstar/wiki/img/GATKstar_flowchart.svg)

## Installation

Clone this repository to your local machine. No additional installation needed.
Please find the information about tool prerequisites at **[Wiki page](https://github.com/gimelbrantlab/GATKstar/wiki)**.

## Citations

Please cite our paper _"Unexpected variability of allelic imbalance estimates from RNA sequencing", Mendelevich A.*, Vinogradova S.*, Gupta S., Mironov A., Synyaev S., Gimelbrant A._, if you used our pipeline in your work.

## License

MIT License







