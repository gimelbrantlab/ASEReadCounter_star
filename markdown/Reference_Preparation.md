# Files needed for the analysis

1. Reference genome

> For example, `GRCm38_68.fa`.
```
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa
```

2. gtf annotation for the reference genome

> For example, `Mus_musculus.GRCm38.68.gtf.gz` for corresponding reference genome version.
```
wget ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
```

3. Either:
* One/Two vcf files for maternal and paternal imbred lines (should be compatible with the reference genome)
* Joint vcf file for multiple species, where two lines are presented (should be compatible with reference genome)
* Individual vcf file or joint individuals vcf file (should be compatible with reference genome)
(see correponding section about reference preparation)
    
> For example, `mgp.v5.merged.snps_all.dbSNP142.vcf.gz` for multiple mice lines.
```
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
```

> *Note: you can do vcf calling by using mutect, varscan, etc. if you heve no pre-existing vcf*


# Input preprocessing before anything:

1. Reference genome fasta should be indexed (`samtools faidx`).

2. vcf-files should be compressed (`bgzip`) and indexed (`samtools tabix`).


# Reference preparation:

One script to rule them all:

```
python3 /full/path/to/ASE/python/prepare_reference_tmp.py --PSEUDOREF True --HETVCF True \
  --pseudoref_dir /full/path/to/dir/for/pseudo/ref/out/ \
  --vcf_dir /full/path/todir/for/vcf/outputs/ \
  --ref /full/path/to/GRCm38_68.fa \
  --name_mat 129S1_SvImJ --name_pat CAST_EiJ \
  --vcf_joint mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
  --gtf /n/scratch2/sv111/ASE/Mus_musculus.GRCm38.68.gtf
```
For help: 
```
python3 /home/am717/ASE/python/prepare_reference_tmp.py --help
```


## Pseudoreference fasta creation:
> `--PSEUDOREF True`

* Input:

|  | Inbred lines | Inbred lines | Individual | Individual | 
| --- | --- | --- | --- | --- |
|  | Joint lines vcf | Separate line(s) vcf | Joint individuals vcf | Separate individual vcf |
| FASTA Reference genome | --ref | --ref | --ref | --ref |
| VCF Variant file(s)[1]    | --vcf_joint | --vcf_mat, --vcf_pat | --vcf_joind | --vcf_ind |
| Name(s)                | --name_mat, --name_pat | --name_mat, --name_pat | --name_ind | --name_ind |
| FASTA Output directory | --pseudo_dir | --pseudo_dir | --pseudo_dir | --pseudo_dir |
| VCF Output directory   | --vcf_dir | --vcf_dir | --vcf_dir | --vcf_dir |

[1] If one or the alleles in case of inbred lines is reference, then everything should be provided as mat or pat only, consistently.

* Output:
  * Pseudoreference genome fastas with own directories.
  * Support vcf or bed files (if needed).


## Heterozygous(parental) VCF creation:
> `--HETVCF True`

* Input:

|  | Inbred lines | Inbred lines | Individual | Individual | 
| --- | --- | --- | --- | --- |
|  | Joint lines vcf | Separate line(s) vcf | Joint individuals vcf | Separate individual vcf |
| FASTA Reference genome | --ref | --ref | --ref | --ref |
| nothing or GTF or BED Selected regions annotation [2] | --gtf or --bed | --gtf or --bed | --gtf or --bed | --gtf or --bed |
| VCF Variant file(s)[1]    | --vcf_joint | --vcf_mat, --vcf_pat | --vcf_joind | --vcf_ind |
| Name(s)                | --name_mat, --name_pat | --name_mat, --name_pat | --name_ind | --name_ind |
| VCF Output directory   | --vcf_dir | --vcf_dir | --vcf_dir | --vcf_dir |

[1] If one or the alleles in case of inbred lines is reference, then everything should be provided as mat or pat only, consistently.
[2] Bed will be considered as main; if gtf provided, automatically considered regions=exons and groups=genes; bed file should have 4 columns and prepared in advance: contig, start position, end position, group ID (no colnames).

* Output:
  * VCF with heterozygous positions, one allele as reference and the second as alternative.
  * Support vcf files (if needed).


# RNA-seq preparation:

## Alignment (STAR) on parental genomes:

* Make sure that your fasta files are ready for the STAR alignment step, each of them should be [indexed with STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf) (`STAR --runMode genomeGenerate`).
```
pseudoRefDirs=/full/path/to/dir/for/pseudo/ref/out/
STAR --runMode genomeGenerate --genomeDir $pseudoRefDirs/129S1_SvImJ/ --genomeFastaFiles $pseudoDir/129S1_SvImJ/129S1_SvImJ_pseudo.fa
STAR --runMode genomeGenerate --genomeDir $pseudoRefDirs/CAST_EiJ/ --genomeFastaFiles $pseudoDir/CAST_EiJ/CAST_EiJ_pseudo.fa
```
* Make sure that you have unziped (`gunzip`) sample fasta/fastq files.

## Creating SNP / Grouped SNPs tables:

Make sure that you have: 
* bed file (with four columns: contig, start and end positions, group ID; no column names) with selected regions, for example, exon regions grouped by genes:
```
awk '$3=="exon" && ($1 ~ /^[1-9]/ || $1 == "X" || $1 == "Y")' /full/path/to/Mus_musculus.GRCm38.68.gtf | cut -f1,4,5,9 | awk -F $'\t' 'BEGIN {OFS = FS} {split($4, a, "gene_id \""); split(a[2], b, "\""); print $1, $2-1, $3, b[1]}' > /full/path/to/output/Mus_musculus.GRCm38.68.EXONS.bed
```
* snp table, which can be obtained from vcf with heterozigous positions created above, as 5 first columns, for example:
```
grep "^#CHROM" /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf | cut -f1-5 > /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt
grep -v "^#" /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf | sort -V | cut -f1-5 >> /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt
```

