#!/usr/bin/env Rscript
# ***********************************************
# Title       : Allelic Imbalance Pipeline
# Description : Script that takes $name_stat_final.txt files as input,
#               merges these files, intersects counts with exons annotaion
#
# Author      : Svetlana Vinogradova
# Date        : 08/24/17
# ***********************************************
#
# run:
#
# Rscript --vanilla counts_to_SNPs_extended2.R -p "/Users/svetlana/Dropbox (Partners HealthCare)/ASE/resultsUpdAugust/Abelson" -r Aleson_test -n Abl1,Abl2 -s GRCm38 -e /Users/svetlana -t exon -m /Users/svetlana/extract_vcf.txt
#
#
# [TODO] : path to lib directory
library("optparse") #, lib.loc="/home/am717/R/x86_64-pc-linux-gnu-library/3.4)"
library("tidyverse") #, lib.loc="/home/am717/R/x86_64-pc-linux-gnu-library/3.4")
#
#
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="path to directory contatining allele counts, output directory would be the same if odir not provided", metavar="character"),
  make_option(c("-n", "--names"), type="character", default=NULL,
              help="names of clones in comma-separated format, for example clone1,clone2,clone3 ", metavar="character"),
  make_option(c("-s", "--suffix_name"), type="character", default=".stat_0.txt",
              help="provide the name suffix for the input files", metavar="character"),
  make_option(c("-r", "--pr_name"), type="character", default="clones",
              help="project name, for example Abelson_clones [default= %default]", metavar="character"),
  make_option(c("-b", "--regions_bed"), type="character", default=NULL,
              help="path to regions of interest bed file (contig, start_pos, end_pos, group_id)", metavar="character"),
  make_option(c("-v", "--snp_tab"), type="character", default=NULL,
              help="path to table with SNPs (from vcf): positions, names, info (first 5 cols)", metavar="character"),
  make_option(c("-o", "--odir"), type="character", default=NULL,
              help="path to output directory", metavar="character")
);
## 0. names and suffixes -- output of allelecout.py
## 1. snp_tab example:
## grep "^#CHROM" $vcfF1 | cut -f1-5 > $snpf1infoexons
## grep -v "^#" $vcfF1 | sort -V | cut -f1-5 > $snpf1infoexons
## [1-based coordinates]
## 2. regions_bed command example (gtf to chromasomal exons bed, with group by gene):
## awk '$3=="exon" && ($1 ~ /^[1-9]/ || $1 == "X" || $1 == "Y")' /n/scratch2/am717/references/hs37d5/Homo_sapiens.GRCh37.63.gtf | cut -f1,4,5,9 | awk -F $'\t' 'BEGIN {OFS = FS} {split($4, a, "gene_id \""); split(a[2], b, "\""); print $1, $2-1, $3, b[1]}' > /n/scratch2/am717/references/hs37d5/Homo_sapiens.GRCh37.63.EXONS.bed
## [0-based coordinates]
#
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
# Check that all arguments are passed
if (is.null(opt$dir) | is.null(opt$names) | is.null(opt$regions_bed) | is.null(opt$snp_tab)){
  print_help(opt_parser)
  stop("Paths and names must be supplied", call.=FALSE)
}
if (is.null(opt$odir)){
  opt$odir <- opt$dir
}
#
# Paths to bed and txt support files (from gtf and vcf):
#
regions_bed <- read_delim(opt$regions_bed, delim="\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, col_types = "ciic")
names(regions_bed) <- c("contig", "start", "end", "groupID")
#
snp_tab <- read_delim(opt$snp_tab, delim="\t", escape_double = FALSE, trim_ws = TRUE, col_types = "ciccc")
colnames(snp_tab) <- c("contig" ,"position", "variantID", "refAllele", "altAllele")
#
print("Contigs in snp file:")
print(unique(snp_tab$contig))
#
# Read input files with SNP allele counts and merge them:
#
names_list <- unlist(strsplit(opt$names, ","))
merged_counts <- snp_tab
for (i in (1:length(names_list))) {
  name <- file.path(opt$dir, paste0(names_list[i], opt$suffix_name))
  rep_tab <- read_delim(name, delim="\t", escape_double = FALSE, trim_ws = TRUE,
                        col_types = cols(contig = col_character()))[ ,c("contig" ,"position", "refCount", "altCount")]
  merged_counts <- merge(merged_counts, rep_tab, by.x=c("contig", "position"), by.y=c("contig", "position"), all.x = TRUE)
  colnames(merged_counts)[1:5] <- c("contig" ,"position", "variantID", "refAllele", "altAllele")
}
merged_counts[is.na(merged_counts)] <- 0
colnames(merged_counts) <- c("contig" ,"position", "variantID", "refAllele", "altAllele",
                             paste0(c("ref", "alt"), "_rep", rep(1:length(names_list), each=2)))
#
print("Contigs after merge:")
print(unique(merged_counts$contig))
#
# Create bed file for bedtools:
#
merged_counts$unique_id <- paste0(merged_counts$contig, "_", merged_counts$position)
# (because vcf is 1-based and bed is 0-based):
merged_counts$start     <- merged_counts$position - 1
merged_counts$end       <- merged_counts$position
#
merged_counts_bed <- merged_counts[,c('contig', 'start', 'end', 'unique_id')]
bed_snp <- file.path(opt$odir, paste0(opt$pr_name, "_merged.v3.1.bed"))
write.table(merged_counts_bed, file=bed_snp, quote = F, row.names = F, col.names = F, sep="\t")
#
# Run bedtools to merge exon data and SNPs:
#
bed_exons_snp <- file.path(opt$odir, paste0(opt$pr_name, "_merged_to_exons.v3.1.bed"))
cmd <- paste("bedtools intersect", "-a", bed_snp, "-b", opt$regions_bed, "-wa", "-wb", ">", bed_exons_snp, sep=" ")
print(paste("CMD::", cmd))
system(cmd)
#
# Merge adding group info (ID) :
#
exons_snp_tab <- read_delim(bed_exons_snp, delim="\t", escape_double = FALSE, trim_ws = TRUE,
                            col_names = FALSE, col_types = "ciicciic")
snp_group_tab <- exons_snp_tab[, c("X4", "X8")]
names(snp_group_tab) <- c("unique_id", "groupID")
merged_counts_group <- merge(merged_counts, snp_group_tab, by="unique_id")
names(merged_counts_group)[1] <- "ID"
merged_counts_group$ID <- paste(merged_counts_group$ID, merged_counts_group$variantID, sep="_")
merged_counts_group <- cbind(merged_counts_group[, !(names(merged_counts_group) %in% c("start", "end", "variantID", "contig",
                                                                                       "position", "refAllele", "altAllele"))],
                             merged_counts_group[, c("contig", "position", "refAllele", "altAllele")])
merged_counts_group <- merged_counts_group[!duplicated(merged_counts_group), ]
mcg_multaffiliation <- merged_counts_group[duplicated(merged_counts_group$ID), ]$ID
merged_counts_group <- merged_counts_group[!(merged_counts_group$ID %in% mcg_multaffiliation), ]
#
# Write SNP allele counts file:
#
file_out_s <- file.path(opt$odir, paste0(opt$pr_name, "_processed_snp.v3.1.txt"))
write.table(merged_counts_group, file=file_out_s, quote = F, row.names = F, sep="\t")
#
# Write gene allele counts file:
#
file_out_g <- file.path(opt$odir, paste0(opt$pr_name, "_processed_gene.v3.1.txt"))
#
merged_counts_group <- merged_counts_group[order(merged_counts_group$contig, merged_counts_group$position), ]
#
gene_snp_num <- aggregate(merged_counts_group[,c(2:(1+(length(names_list)*2)))], by=list(merged_counts_group$groupID), length)[, c(1:2)]
names(gene_snp_num) <- c("ID", "n_snp")
gene_snp_pos <- aggregate(merged_counts_group[,c("contig", "position")], by=list(merged_counts_group$groupID), paste, collapse=',')
names(gene_snp_pos) <- c("ID", "contig", "aggr_pos")
gene_snp_pos$contig <- sapply(gene_snp_pos$contig, function(x){strsplit(x, ',')[[1]][1]})
#
gene_snp_sum <- aggregate(merged_counts_group[, c(2:(1+(length(names_list)*2)))],
                          by=list(merged_counts_group$groupID), sum)
names(gene_snp_sum)[1] <- "ID"
#
gene_snp_aggr <- aggregate(merged_counts_group[,c(2:(1+(length(names_list)*2)))], by=list(merged_counts_group$groupID), paste, collapse=",")
names(gene_snp_aggr)[1] <- "ID"
names(gene_snp_aggr)[2:ncol(gene_snp_aggr)] <- paste0("aggr_", names(gene_snp_aggr)[2:ncol(gene_snp_aggr)])
#
# MERGE:
#
merge_gene_tab <- merge(gene_snp_sum, gene_snp_num, by="ID")
merge_gene_tab <- merge(merge_gene_tab, gene_snp_pos, by="ID")
merge_gene_tab <- merge(merge_gene_tab, gene_snp_aggr, by="ID")
merge_gene_tab <- merge_gene_tab[order(merge_gene_tab$contig,
                                       as.numeric(sapply(merge_gene_tab$aggr_pos, function(x){strsplit(x, ',')[[1]][1]}))), ]
#
write.table(merge_gene_tab, file=file_out_g, sep = "\t", quote = F, row.names = F)
#
#
#
#
#
#
names(merge_gene_tab)[1] = "ensembl_gene_id"
names(merge_gene_tab)[names(merge_gene_tab)=="contig"] = "chr"
write.table(merge_gene_tab, file=file.path(opt$odir, paste0(opt$pr_name, "_processed_gene_extended2.txt")),
            sep = "\t", quote = F, row.names = F)
#


































