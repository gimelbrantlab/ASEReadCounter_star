
#########################################################
# 							#
# 		ASE Replicates Project			#
#         https://github.com/gimelbrantlab/ASE		#
#							#
#   Authors: Mendelevich Asya, Svetlana Vinogradova	#
#							#
#########################################################
#
# DESCRIPTION:
#   Function for full reference preprocessing:
#     - Creating pseudogenomes from reference genomes and vcf file(s)
#     - F1-cross vcf files
#     - Gene-Transcript-Exon Annotations
# 
# NOTE: the order of chromosomes in each file should be the same!
#
# DEPEND:
# bash(v4.2.46)
# gcc(v6.2.0)
# python(v3.6.0)
# java(v1.8.0_112)
# gatk(v4.0.0.0)
# bcftools(v1.3.1)
# htslib(v1.3.2)
# cufflinks(v2.2.1)
# vcftools(v0.1.17)
#
# Please, use --help|-h for help. 
# For more information, please, visit:
#       https://github.com/gimelbrantlab/GATKstar
#

import argparse
import os
import tempfile
import subprocess
import gzip
import re


def GATK_SelectVariants(r, v, o, g=None, n=None, b=False):
    '''
    The GATK command for VCF processings, selects SNPs from all variants
    (optionally, restricts on exons only, or takes only biallelic variants)
    Input:  
    r -- genome.fa (should be indexed (samtools faidx) and have dictionary file (Picard CreateSequenceDictionary))
    v -- vcf file (with one column and with only biallelic variants for particular organism (not a set) remains)
    o -- ofile to place output 
    g -- (optional) gtf file for exon positions annotation
    n -- (optional) name of the column to chop from mixed vcf 
    b -- (optional) restriction to biallelic variants flag
    Output: 
    returns nothing
    creates vcf file (name defined via o option) with selected variants
    '''

    cmd = "gatk SelectVariants -select-type SNP "
    flags = {'-R':r, '-V':v, '-O':o}
    if (g):
        # tmp exons bed file:
        #exon_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
        exon_bed = os.path.join(os.path.dirname(o), os.path.basename(g) + "exons.bed")
        cmd_exon = " ".join(["grep -w 'exon'", g, "| grep '^[0-9XY]' | awk 'BEGIN{FS=OFS=", '"\t"', "}; {print $1,$4-1,$5}' >", exon_bed])
        print(cmd_exon)
        subprocess.check_output(cmd_exon, shell=True)
        flags['-L'] = exon_bed
    if (n):
        flags['-sn'] = n
    if (b):
        flags['--restrict-alleles-to'] = "BIALLELIC"

    flags_str = ''
    for f in flags:
        if isinstance(flags[f], str):
            flags_str += f + ' ' + flags[f] + ' '
        else:
            for item in flags[f]:
                flags_str += f + ' ' + item + ' '
    cmd = ' '.join([cmd, flags_str])
    print(cmd)
    subprocess.check_output(cmd, shell=True)

    # clear up!:
    #if (g):
    #    os.remove(exon_bed.name)
    return

def SelectBiallelicSNP_VCF(v, o, name):
    '''
    GATK biallelic variant selection does not work properly.
    That's a replacement. 
    Input:
    v    -- path to SNP vcf
    o    -- path to output biallelic SNP vcf
    name -- name of the column
    Output: 
    vcf file with biallelic (REF!=ALT) sites. 
    '''
    vcf_stream = open(v, 'r') 
    out_stream = open(o, 'w')

    # Read header:
    row = vcf_stream.readline()
    while (row.startswith("##")):
        out_stream.write(row)
        row = vcf_stream.readline()

    # Column Names: 
    colnames = row.replace('#','').strip().split()
    name_col = colnames.index(name)
    format_col = colnames.index("FORMAT")
    ref_col = colnames.index("REF")
    alt_col = colnames.index("ALT")
    out_stream.write(row)

    # Row by row:
    for row in vcf_stream:
        row = row.strip().split()
        gt_index = row[format_col].split(":").index("GT")
        gt_name_list = row[name_col].split(":")
        gt_name = row[name_col].split(":")[gt_index]

        if (len(gt_name)==1 and gt_ind!='.'):

            if (gt_name != '0'):
                alt_allele = row[alt_col].split(',')[int(gt_name)-1]
                gt_name = '1|1'

                gt_name_list[gt_index] = gt_name
                row[name_col] = ":".join(gt_name_list)
                row[alt_col] = alt_allele
                row_out = '\t'.join(row)
                out_stream.write(row_out + '\n')

        elif (gt_name[0]!='.' and gt_name[2]!='.'):

            if (gt_name[0] != '0' and gt_name[0] == gt_name[2]):
                alt_allele = row[alt_col].split(',')[int(gt_name[0])-1]
                gt_name = '1|1'

                gt_name_list[gt_index] = gt_name
                row[name_col] = ":".join(gt_name_list)
                row[alt_col] = alt_allele
                row_out = '\t'.join(row)
                out_stream.write(row_out + '\n')

    vcf_stream.close()
    out_stream.close()
    return


def ParentalSeparation_VCF(v, o_1, o_2, ind_name):
    '''
    Takes individual fazed vcf and splits it in two ref-alt1/alt2 files.
    Input:
    v  -- path to individual vcf
    o1 -- path to output file for first allele vcf
    o2 -- path to output file for second allele vcf
    ind_name -- name of individual, should coinside with column name in vcf
    Output:
    Two vcf files for two alleles of the organism, reference is reference, haplotype 1|1, for each allele
    '''
    vcf_stream = open(v, 'r')
    out1_stream = open(o_1, 'w')
    out2_stream = open(o_2, 'w')

    # Read header:
    row = vcf_stream.readline()
    while (row.startswith("##")):
        out1_stream.write(row)
        out2_stream.write(row)
        row = vcf_stream.readline()

    # Column Names: 
    colnames = row.replace('#','').strip().split()
    ind_col = colnames.index(ind_name)
    format_col = colnames.index("FORMAT")
    ref_col = colnames.index("REF")
    alt_col = colnames.index("ALT")

    colnames[ind_col] = ind_name + ".mat"
    out1_stream.write('#' + '\t'.join(colnames) + '\n')
    colnames[ind_col] = ind_name + ".pat"
    out2_stream.write('#' + '\t'.join(colnames) + '\n')

    # Row by row:
    for row in vcf_stream:
        row = row.strip().split()
        gt_index = row[format_col].split(":").index("GT")
        gt_ind_list = row[ind_col].split(":")
        gt_ind = row[ind_col].split(":")[gt_index]

        if (len(gt_ind)==1 and gt_ind!='.'):
            
            if (gt_ind == '0'):
                gt_ind_12 = '0|0'
                alt_allele_12 = row[ref_col]
            else:
                gt_ind_12 = '1|1'
                alt_allele_12 = row[alt_col].split(',')[int(gt_ind)-1]

            gt_ind_list[gt_index] = gt_ind_12
            row[ind_col] = ":".join(gt_ind_list)
            row[alt_col] = alt_allele_12
            cols_out12 = '\t'.join(row)
            out1_stream.write(cols_out12 + '\n')
            out2_stream.write(cols_out12 + '\n')

        elif (gt_ind[0]!='.' and gt_ind[2]!='.'):

            if (gt_ind[0] == '0'):
                gt_ind_1 = '0|0'
                alt_allele_1 = row[ref_col]
            else:
                gt_ind_1 = '1|1'
                alt_allele_1 = row[alt_col].split(',')[int(gt_ind[0])-1]

            if (gt_ind[2] == '0'):
                gt_ind_2 = '0|0'
                alt_allele_2 = row[ref_col]
            else:
                gt_ind_2 = '1|1'
                alt_allele_2 = row[alt_col].split(',')[int(gt_ind[2])-1]

            gt_ind_list[gt_index] = gt_ind_1
            row[ind_col] = ":".join(gt_ind_list)
            row[alt_col] = alt_allele_1
            cols_out1 = '\t'.join(row)
            out1_stream.write(cols_out1 + '\n')

            gt_ind_list[gt_index] = gt_ind_2
            row[ind_col] = ":".join(gt_ind_list)
            row[alt_col] = alt_allele_2
            cols_out2 = '\t'.join(row)
            out2_stream.write(cols_out2 + '\n')

        ### IF THERE ARE OTHER FIELDS IN COLUMN TO BE CHANGED? ###

    vcf_stream.close()
    out1_stream.close()
    out2_stream.close()
    return


def gzip_tabix_VCF(vcf):
    '''bgzip+tabix'''
    cmd_bgzip = " ".join(["bgzip -c", vcf, ">", vcf + '.gz'])
    print(cmd_bgzip)
    subprocess.check_output(cmd_bgzip, shell=True)
    cmd_tabix = " ".join(["tabix -p vcf", vcf + '.gz']),
    print(cmd_tabix)
    subprocess.check_output(cmd_tabix, shell=True)
    return

def vcftools_consensus(ref_fa, vcf, pseudo_fa):
    '''creates pseudogenom inserting corresponding SNPs (vcftools vcf-consensus)'''
    # it is VERY SLOW, maybe replace with something else?
    #cmd = " ".join(["<", ref_fa, "vcf-consensus", vcf, ">", pseudo_fa])
    cmd = " ".join(["cat", ref_fa, "| vcf-consensus", vcf, ">", pseudo_fa])
    print(cmd)
    subprocess.check_output(cmd, shell=True)
    return


def SepToPair_VCF(ref, vcf_mat, vcf_pat, ofile, name_mat, name_pat):
    unfiltered = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
    cmd_merge = " ".join(["bcftools merge", vcf_mat, vcf_pat, "-m both -o", unfiltered.name])
    print(cmd_merge)
    subprocess.check_output(cmd_merge, shell=True)
    GATK_SelectVariants(r=ref, v=unfiltered.name, o=ofile)
    os.remove(unfiltered.name)
    return

def PairToF1_VCF(vcf_pair, vcf_f1, name_mat, name_pat):
    print(vcf_pair + "  -->  " + vcf_f1)
    vcf_stream = open(vcf_pair, 'r')
    out_stream = open(vcf_f1, 'w')

    # header:
    ### REWRITE ACCURATE header filtering ###
    row = vcf_stream.readline()
    while (row.startswith("##")):
        out_stream.write(row)
        row = vcf_stream.readline()

    colnames = row.replace('#','').strip().split()
    mat_col = colnames.index(name_mat)
    pat_col = colnames.index(name_pat)
    format_col = colnames.index("FORMAT")
    ref_col = colnames.index("REF")
    alt_col = colnames.index("ALT")

    colnames[pat_col] = "F1"
    colnames.pop(mat_col)
    out_stream.write('#' + '\t'.join(colnames) + '\n')

    # body:
    for row in vcf_stream:
        row = row.strip().split('\t')
        gt_index = row[format_col].split(":").index("GT")
        gt_mat_list = row[mat_col].split(":")
        gt_mat = row[mat_col].split(":")[gt_index]
        gt_pat_list = row[pat_col].split(":")
        gt_pat = row[pat_col].split(":")[gt_index]

        if (len(gt_mat)==1):
            gt_mat = gt_mat + '|' + gt_mat
        if (len(gt_pat)==1):
            gt_pat = gt_pat + '|' + gt_pat

        if (gt_mat[0] != '.' and gt_pat[0] != '.' and gt_mat[0]==gt_mat[2] and gt_pat[0]==gt_pat[2] and gt_mat[0]!=gt_pat[0]):
            if (gt_mat[0] == '0'):
                ref_allele = row[ref_col]
            else:
                ref_allele = row[alt_col].split(',')[int(gt_mat[0])-1]
            if (gt_pat[0] == '0'):
                alt_allele = row[ref_col]
            else:
                alt_allele = row[alt_col].split(',')[int(gt_pat[0])-1]

            gt_pat_list[gt_index] = "0|1"
            row[pat_col] = ":".join(gt_pat_list)
            row[ref_col] = ref_allele
            row[alt_col] = alt_allele
            row.pop(mat_col)

            ### IF THERE ARE OTHER FIELDS IN COLUMN TO BE CHANGED? ###

            row_out = '\t'.join(row)
            out_stream.write(row_out + '\n')

    vcf_stream.close()
    out_stream.close()
    return

def SingletonToF1_VCF(vcf_single, vcf_f1, name_alt):
    print(vcf_single + "  -->  " + vcf_f1)
    vcf_stream = open(vcf_single, 'r')
    out_stream = open(vcf_f1, 'w')

    # header:
    ### REWRITE ACCURATE header filtering ###
    row = vcf_stream.readline()
    while (row.startswith("##")):
        out_stream.write(row)
        row = vcf_stream.readline()

    colnames = row.replace('#','').strip().split()
    name_col = colnames.index(name_alt)
    format_col = colnames.index("FORMAT")
    ref_col = colnames.index("REF")
    alt_col = colnames.index("ALT")

    colnames[name_col] = "F1"
    out_stream.write('#' + '\t'.join(colnames) + '\n')

    # body:
    for row in vcf_stream:
        row = row.strip().split('\t')
        gt_index = row[format_col].split(":").index("GT")
        gt_name_list = row[name_col].split(":")
        gt_name = row[name_col].split(":")[gt_index]

        if (len(gt_name)==1):
            gt_name = gt_name + '|' + gt_name

        if (gt_name[0]==gt_name[2] and gt_name[0]!='0' and gt_name[0]!='.'):
            row[alt_col] = row[alt_col].split(',')[int(gt_name[0])-1]

            gt_name_list[gt_index] = "0|1"
            row[name_col] = ":".join(gt_name_list)

            ### IF THERE ARE OTHER FIELDS IN COLUMN TO BE CHANGED? ###

            row_out = '\t'.join(row)
            out_stream.write(row_out + '\n')

    vcf_stream.close()
    out_stream.close()
    return


def IndToF1_VCF(vcf_ind, vcf_f1, name_ind):
    print(vcf_ind + "  -->  " + vcf_f1)
    vcf_stream = open(vcf_ind, 'r')
    out_stream = open(vcf_f1, 'w')

    row = vcf_stream.readline()
    while (row.startswith("##")):
        out_stream.write(row)
        row = vcf_stream.readline()

    colnames = row.replace('#','').strip().split()
    ind_col = colnames.index(name_ind)
    format_col = colnames.index("FORMAT")
    ref_col = colnames.index("REF")
    alt_col = colnames.index("ALT")

    colnames[ind_col] = "F1"
    out_stream.write('#' + '\t'.join(colnames) + '\n')

    # body:
    for row in vcf_stream:
        row = row.strip().split('\t')
        gt_index = row[format_col].split(":").index("GT")
        gt_name_list = row[ind_col].split(":")
        gt_name = row[ind_col].split(":")[gt_index]

        if (len(gt_name)==3 and gt_name[0]!=gt_name[2] and gt_name[0]!='.' and gt_name[2]!='.'):
            if (gt_name[0] == '0'):
                ref_var = row[ref_col]
            else:
                ref_var = row[alt_col].split(',')[int(gt_name[0])-1]

            if (gt_name[2] == '0'):
                alt_var = row[ref_col]
            else:
                alt_var = row[alt_col].split(',')[int(gt_name[2])-1]

            row[ref_col] = ref_var
            row[alt_col] = alt_var
            gt_name_list[gt_index] = "0|1"
            row[ind_col] = ":".join(gt_name_list)

            ### IF THERE ARE OTHER FIELDS IN COLUMN TO BE CHANGED? ###

            row_out = '\t'.join(row)
            out_stream.write(row_out + '\n')

    vcf_stream.close()
    out_stream.close()
    return



def main():
    # ARGUMENTS HANDLING: 
    parser = argparse.ArgumentParser()
    parser.add_argument("--PSEUDOREF", required=False, default="False", metavar="[True/False]", help="Pseudogenomes creation is needed? Default: False")
    parser.add_argument("--HETVCF", required=False, default="False", metavar="[True/False]", help="Het vcf creation is needed? Default: False")
    parser.add_argument("--pseudoref_dir", required=False, metavar="[/path/to/dir]", help="Path to directory with subdirectories --name_mat and --name_pat, that contain pseudogenome fasta files.")
    parser.add_argument("--vcf_dir", required=False, metavar="[/path/to/dir]", help="Path to directory with vcf files; required if --HETVCF True, or allelevcfs not provided")
    parser.add_argument("--ref", required=True, metavar="[/path/to/file]", help="Path to reference fasta file.")
    parser.add_argument("--gtf", required=False, metavar="[/path/to/file]", help="Path to reference gtf file, if provided will create *.exon.vcf.")
    parser.add_argument("--bed", required=False, metavar="[/path/to/file]", help="Path to bed file with regions for selection (3 first columns (no colnames!): chrom, chromStart, chromEnd), if provided will create *.selected_regions.vcf")
    parser.add_argument("--name_ind", required=False, help="Name of individuum (required with --vcf_joind and --vcf_ind: name should coincide with name in vcf; please avoid giving 'ref' ar 'alt' names, they have special meanings)")
    parser.add_argument("--name_mat", required=False, help="Name of maternal line/sample (required with --vcf_joint if aat is needed, and --vcf_mat: names should coincide with names in vcf; please avoid giving 'ref' ar 'alt' names, they have special meanings)")
    parser.add_argument("--name_pat", required=False, help="Name of paternal line/sample (required with --vcf_joint if pat is needed, and --vcf_pat: names should coincide with names in vcf; please avoid giving 'ref' ar 'alt' names, they have special meanings)")
    parser.add_argument("--vcf_mat", required=False, metavar="[/path/to/file]", help="Path to maternal vcf file (if F1 cross, eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_pat", required=False, metavar="[/path/to/file]", help="Path to paternal vcf file (if F1 cross, eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_joint", required=False, metavar="[/path/to/file]", help="Path to joint vcf file for lines (if F1 cross, eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_ind", required=False, metavar="[/path/to/file]", help="Path to individuum vcf file (if not F1 cross, eigther --vcf_joind or --vcf_ind required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_joind", required=False, metavar="[/path/to/file]", help="Path to joint vcf file for individuums (if not F1 cross, eigther --vcf_joind or --vcf_ind required; separate vcf will dominate if both provided)") 

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # TEST IF EVERYTHING NECESSARY IS PRESENT and SET NAMES and A MODE |
    # ------------------------------------------------------------------
    # GENERAL PARAMs:

    if (args.ref is None):
         msg = "Required parameter --ref is missing."
         raise argparse.ArgumentTypeError(msg)
    if (args.PSEUDOREF is True and args.pseudoref_dir is None):
         msg = "Required parameter --pseudoref_dir is missing."
         raise argparse.ArgumentTypeError(msg)
    # if (args.HETVCF is True and args.vcf_dir is None):
    #      msg = "Required parameter --vcf_dir is missing."
    #      raise argparse.ArgumentTypeError(msg)
    if (args.vcf_dir is None):
         msg = "Required parameter --vcf_dir is missing."

    # CASES:
    # SEPARATE ALLELE VCFs:
    if (args.vcf_mat is not None or args.vcf_pat is not None):
        if ((args.vcf_mat is not None and args.name_mat is not None) and (args.vcf_pat is not None and args.name_pat is not None)):
            input_case = "two_alleles"
            name_mat   = args.name_mat
            name_pat   = args.name_pat
        elif (args.vcf_mat is not None and args.name_mat is not None):
            input_case = "one_allele"
            vcf_alt    = args.vcf_mat
            name_alt   = args.name_mat
        elif (args.vcf_pat is not None and args.name_pat is not None):
            input_case = "one_allele"
            vcf_alt    = args.vcf_pat
            name_alt   = args.name_pat
        else: 
            msg = "Data required: --vcf_xxx should go in a pair with --name_xxx."
            raise argparse.ArgumentTypeError(msg)
    # JOINT VCF FILE WITH ALLELES as columns:
    elif (args.vcf_joint is not None):
        if (args.name_mat is not None and args.name_pat is not None):
            input_case = "joint_two_alleles"
            name_mat   = args.name_mat
            name_pat   = args.name_pat
        elif (args.name_mat is not None):
            input_case = "joint_one_allele"
            name_alt   = args.name_mat
        elif (args.name_pat is not None):
            input_case = "joint_one_allele"
            name_alt   = args.name_pat
        else:
            msg = "Data required: at least either --name_mat or --name_pat is needed."
            raise argparse.ArgumentTypeError(msg)
    # INDIVIDUAL or JOINT INDUVUDUAL VCF (not a cross):
    elif (args.vcf_ind is not None or args.vcf_joind is not None):
        if (args.name_ind is not None):
            name_ind   = args.name_ind
            name_mat   = str(args.name_ind) + "_mat"
            name_pat   = str(args.name_ind) + "_pat"
            if (args.vcf_ind is not None):
                input_case = "individ"
            elif (args.vcf_joind is not None):
                input_case = "joint_individs"
        else:
            msg = "Data required: parameter --name_ind is missing."
            raise argparse.ArgumentTypeError(msg)
    # SMTH INCORRECT:
    else:
        msg = "Check the correctness of required parameters for your case. See help."
        parser.print_help()
        raise argparse.ArgumentTypeError(msg)

    # ------------------------------------------------------------------
    # PSEUDOREFERENCE CREATION: if PSEUDOREF set to be True            |
    # ------------------------------------------------------------------
 
    if (args.PSEUDOREF=="True"):

        #cmd_mkdir_vcf = "mkdir -p " + args.vcf_dir
        #subprocess.check_output(cmd_mkdir_vcf, shell=True)
        os.makedirs(args.vcf_dir, exist_ok=True)

        if (input_case == "two_alleles" or input_case == "joint_two_alleles" or input_case == "individ" or input_case == "joint_individs"):
            # Output vcfs:
            sep_vcf_mat = os.path.join(args.vcf_dir, name_mat + ".SNP.biallelic.vcf")
            sep_vcf_pat = os.path.join(args.vcf_dir, name_pat + ".SNP.biallelic.vcf")
            # Creation of directories for pseudoref:
            pseudo_dir_mat = os.path.join(args.pseudoref_dir, name_mat)
            pseudo_dir_pat = os.path.join(args.pseudoref_dir, name_pat)
            os.makedirs(pseudo_dir_mat, exist_ok=True)
            os.makedirs(pseudo_dir_pat, exist_ok=True)
            #cmd_mkdir_mat = "mkdir -p " + pseudo_dir_mat
            #subprocess.check_output(cmd_mkdir_mat, shell=True)
            #cmd_mkdir_pat = "mkdir -p " + pseudo_dir_pat
            #subprocess.check_output(cmd_mkdir_pat, shell=True)
            # All cases:
            if (input_case == "two_alleles"):
                snp_mat_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
                snp_pat_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")

                GATK_SelectVariants(r=args.ref, v=args.vcf_pat, o=snp_pat_vcf.name, b=False)
                GATK_SelectVariants(r=args.ref, v=args.vcf_mat, o=snp_mat_vcf.name, b=False)
                SelectBiallelicSNP_VCF(snp_mat_vcf.name, sep_vcf_mat, name_mat)
                SelectBiallelicSNP_VCF(snp_pat_vcf.name, sep_vcf_pat, name_pat)
                
                os.remove(snp_mat_vcf.name)
                os.remove(snp_pat_vcf.name)

            elif (input_case == "joint_two_alleles"):
                snp_mat_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
                snp_pat_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
                
                GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=snp_pat_vcf.name, n=name_pat, b=False)
                GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=snp_mat_vcf.name, n=name_mat, b=False)
                SelectBiallelicSNP_VCF(snp_mat_vcf.name, sep_vcf_mat, name_mat)
                SelectBiallelicSNP_VCF(snp_pat_vcf.name, sep_vcf_pat, name_pat)

                os.remove(snp_mat_vcf.name)
                os.remove(snp_pat_vcf.name)

            elif (input_case == "individ" or input_case == "joint_individs"):
                sep_vcf_ind = os.path.join(args.vcf_dir, name_ind + ".individual.SNP.vcf")
                if (input_case == "individ"):
                    GATK_SelectVariants(r=args.ref, v=args.vcf_ind, o=sep_vcf_ind, b=False)
                elif (input_case == "joint_individs"):
                    GATK_SelectVariants(r=args.ref, v=args.vcf_joind, o=sep_vcf_ind, n=name_ind, b=False)
                ParentalSeparation_VCF(v=sep_vcf_ind, o_1=sep_vcf_mat, o_2=sep_vcf_pat, ind_name=name_ind)

            # Indexing vcfs:
            gzip_tabix_VCF(sep_vcf_mat)
            gzip_tabix_VCF(sep_vcf_pat)
            # Pseudoreference:
            vcftools_consensus(args.ref, sep_vcf_mat + '.gz', os.path.join(pseudo_dir_mat, name_mat + "_pseudo.fa"))
            vcftools_consensus(args.ref, sep_vcf_pat + '.gz', os.path.join(pseudo_dir_pat, name_pat + "_pseudo.fa"))

        elif (input_case == "one_allele" or input_case == "joint_one_allele"):
            # Alt vcfs:
            sep_vcf_alt = os.path.join(args.vcf_dir, name_alt + ".SNP.biallelic.vcf")
            # Creation of directories for pseudoref:
            pseudo_dir_alt = os.path.join(args.pseudoref_dir, name_alt)
            os.makedirs(pseudo_dir_alt, exist_ok=True)
            #cmd_mkdir_alt = "mkdir -p " + pseudo_dir_pat
            #subprocess.check_output(cmd_mkdir_alt, shell=True)
            # All cases:
            if(input_case == "one_allele"):
                snp_alt_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
                
                GATK_SelectVariants(r=args.ref, v=args.vcf_mat, o=snp_alt_vcf.name, b=False)
                SelectBiallelicSNP_VCF(snp_alt_vcf.name, sep_vcf_alt, name_alt)

                os.remove(snp_alt_vcf.name)
            
            elif(input_case == "joint_one_allele"):
                snp_alt_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")

                GATK_SelectVariants(r=args.ref, v=args.vcf_mat, n=name_alt, o=snp_alt_vcf.name, b=False)
                SelectBiallelicSNP_VCF(snp_alt_vcf.name, sep_vcf_alt, name_alt)

                os.remove(snp_alt_vcf.name)
            # Indexing vcfs:
            gzip_tabix_VCF(sep_vcf_alt)
            # Pseudoreference:
            vcftools_consensus(args.ref, sep_vcf_alt + '.gz', os.path.join(pseudo_dir_alt, name_alt + "_pseudo.fa"))

        else:
            msg = "Something went wrong in the input data parsing."
            raise argparse.ArgumentTypeError(msg)

    # ------------------------------------------------------------------
    # HETEROZYGOUS SNP VCF CREATION: if HETVCF set to be True          |
    # ------------------------------------------------------------------

    if (args.HETVCF=="True"):

        os.makedirs(args.vcf_dir, exist_ok=True)

        # F1 VCF:

        if (input_case == "two_alleles" or input_case == "joint_two_alleles"):
            vcf_het = os.path.join(args.vcf_dir, "_".join(["Het_Allelic", name_mat, name_pat])+'.vcf')

            # .. to pair (snp):
            pair_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
            if (input_case == "two_alleles"):
                SepToPair_VCF(args.ref, args.vcf_mat, args.vcf_pat, pair_vcf.name, name_mat, name_pat)
            elif (input_case == "joint_two_alleles"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=pair_vcf.name, n=[name_mat, name_pat])

            # pair (snp) to f1:
            PairToF1_VCF(pair_vcf.name, vcf_het, name_mat, name_pat)

            os.remove(pair_vcf.name)

        elif (input_case == "one_allele" or input_case == "joint_one_allele"):
            vcf_het = os.path.join(args.vcf_dir, "_".join(["Het_Allelic", "ref", name_alt])+'.vcf')

            # .. to singleton (snp):
            singleton_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
            if (input_case == "one_allele"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_alt, o=singleton_vcf.name)
            elif (input_case == "joint_one_allele"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=singleton_vcf.name, n=name_alt)          

            # singleton (snp) to f1:
            SingletonToF1_VCF(singleton_vcf.name, vcf_het, name_alt)
            
            os.remove(singleton_vcf.name)

        elif (input_case == "individ" or input_case == "joint_individs"):
            vcf_het = os.path.join(args.vcf_dir, "_".join(["Het_Allelic", name_ind])+'.vcf')

            # .. to sep ind vcf (snp):
            sep_vcf_ind = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
            if (input_case == "individ"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_ind, o=sep_vcf_ind.name)
            elif (input_case == "joint_individs"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_joind, o=sep_vcf_ind.name, n=name_ind)

            # sep ind vcf (snp) to f1:
            IndToF1_VCF(sep_vcf_ind.name, vcf_het, name_ind)

            os.remove(sep_vcf_ind.name)

        gzip_tabix_VCF(vcf_het)

        # HET VCF restrictions on the regions:
        if (args.bed is not None or args.gtf is not None):
            if (args.bed is not None):
                vcf_het_region = vcf_het.replace(".vcf", ".selected_regions.vcf")
                
                region_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")

                cmd_bed_head = "sed '1ichrom\tchromStart\tchromEnd' " + args.bed + " > " + region_bed 
                print(cmd_bed_head)
                subprocess.check_output(cmd_bed_head, shell=True)

            elif (args.gtf is not None):
                vcf_het_region = vcf_het.replace(".vcf", ".exons.vcf")

                exon_bed = os.path.join(os.path.dirname(vcf_het_exon), os.path.basename(args.gtf) + "exons.bed")
                cmd_exon_bed = " ".join(["grep -w 'exon'", args.gtf, "| grep '^[0-9XY]' | awk 'BEGIN{FS=OFS=", '"\t"', "}; {print $1,$4-1,$5}' >", exon_bed])
                print(cmd_exon_bed)
                subprocess.check_output(cmd_exon_bed, shell=True)
                cmd_bed_head = "sed -i '1ichrom\tchromStart\tchromEnd' " + exon_bed
                print(cmd_bed_head)
                subprocess.check_output(cmd_bed_head, shell=True)
   
                region_bed = exon_bed
 
            cmd_getregion_vcf = ' '.join(["vcftools", "--vcf", vcf_het, "--bed", region_bed, "--recode", "--recode-INFO-all", "--stdout", '>', vcf_het_region])
            print(cmd_getregion_vcf)
            subprocess.check_output(cmd_getregion_vcf, shell=True)
 
            gzip_tabix_VCF(vcf_het_region)

            if (args.bed is not None):
                os.remove(region_bed)




if __name__ == "__main__":
    main()

