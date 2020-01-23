# GATK* pipeline for Allelic Counting
([GATK](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)-like) data preprocessing for Allele-specific expression analysis on RNA-seq data.

This pipeline was used as the main procedure for creating tables of gene allelic counts for the further AI analysis with [QCumber](https://github.com/gimelbrantlab/QCumber) in paper

>_"Unexpected variability of allelic imbalance estimates from RNA sequencing", Mendelevich A.*, Vinogradova S.*, Gupta S., Mironov A., Synyaev S., Gimelbrant A._

The pipeline can be divided into two prime parts:

* **Reference preparation**
  constructs individual paternal and maternal genomes, creates heterozigous VCF.

* **Creation of tables with Allelic Counts**
  maps the reads from RNA-seq experiments to these genomes and counts the number of reads which map to either the reference or alternate allele at each heterozygous SNP, estimates allelic imbalance for individual genes summarizing information from SNPs.

Please find **manuals / worked examples** at **[Wiki page](https://github.com/gimelbrantlab/GATKstar/wiki)**.

The resulting allelic counts tables can be used for allelic imbalance analysis via [QCumber](https://github.com/gimelbrantlab/QCumber).


![scheme]()



## Credits
TODO: Write credits

## License
TODO: Write license








