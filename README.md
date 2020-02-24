# ASEReadCounter* - preprocessing sequencing data for allele-specific analysis
This pipeline goes from RNA-seq (or similar) data to a table of total allelic counts per gene (or other genomic interval). That table serves as input for the further analysis of allelic imbalance with [Qllelic](https://github.com/gimelbrantlab/Qllelic). 

This is a re-implementation of the `ASEReadCounter` tool from [GATK](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6), based on [allelecounter](https://github.com/secastel/allelecounter) scripts by S.Castel.

The pipeline consists of two main parts:

0. **Reference preparation**

     construct individual "paternal" and "maternal" genome references, create heterozygous VCF.

1. **Creation of tables with allelic counts**

  * map sequencing reads to references (using STAR aligner; see [complete list of dependencies](https://github.com/gimelbrantlab/ASEReadCounter_star/wiki/Prerequisites:-tools-and-packages))
  * perform random sampling of the mapped reads to defined depth (key step for overdispersion analysis in [Qllelic](https://github.com/gimelbrantlab/Qllelic))
  * count the number of reads mapping to the reference or alternate allele at each heterozygous SNP, and collate the counts for genome intervals (e.g., genes or other features).

Please find **manuals / worked examples** at **[Wiki page](https://github.com/gimelbrantlab/ASEReadCounter_star/wiki)** of this repository.


![pic](https://github.com/gimelbrantlab/ASEReadCounter_star/blob/master/ASEReadCounter_star_flowchart_short.svg)


## Installation

Clone this repository to your local machine. No additional installation needed.
Please find the information about tool prerequisites at **[Wiki page](https://github.com/gimelbrantlab/ASEReadCounter_star/wiki)**.

## Citations

Please cite _"Unexpected variability of allelic imbalance estimates from RNA sequencing", Mendelevich A.*, Vinogradova S.*, Gupta S., Mironov A., Sunyaev S., Gimelbrant A._, if you used our pipeline in your work.

## Reporting bugs

Please report bugs to the Github [issues](https://github.com/gimelbrantlab/ASEReadCounter_star/issues) page.

## License

[GNU General Public License v3.0](https://github.com/gimelbrantlab/ASEReadCounter_star/blob/master/LICENSE)








