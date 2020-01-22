##
## NOTE: 
## -- SAM-files should be already sorted by read names
## -- Read names must be non-empty
## USAGE: alleleseq_merge_stream.py --pat_sam [path to file] --mat_sam [path to file] --o [path to file] --paired [0|1]
##

#import cProfile
import sys
import pandas
import argparse
import random
import subprocess
import time
start_time = time.time()

# special samtools ordering:
def isdigit(c):
    return c.isdigit()
def uporotiycmp(firstrb, secstrb):
    firstr = firstrb + '\x00' 
    secstr = secstrb + '\x00'
    firind , secind = 0 , 0
    firlen = len(firstr)
    seclen = len(secstr)
    while (firind < firlen) and (secind < seclen):
        if (isdigit(firstr[firind]) and isdigit(secstr[secind])):
            while (firstr[firind] == '0'):
                firind += 1
            while (secstr[secind] == '0'):
                secind += 1
            while isdigit(firstr[firind]) and isdigit(secstr[secind]) and (firstr[firind] == secstr[secind]):
                firind += 1
                secind += 1
            if (isdigit(firstr[firind]) and isdigit(secstr[secind])):
                i = 0
                while (isdigit(firstr[firind+i]) and isdigit(secstr[secind+i])):
                    i += 1
                if isdigit(firstr[firind+i]): return 1
                elif isdigit(secstr[secind+i]): return -1
                else: return ord(firstr[firind]) - ord(secstr[secind])
            elif (isdigit(firstr[firind])): return 1
            elif (isdigit(secstr[secind])): return -1
            elif (firind < secind): return 1
            elif (firind > secind): return -1
        else:
            if (firstr[firind] != secstr[secind]):
               return ord(firstr[firind]) - ord(secstr[secind])
            firind += 1
            secind += 1
    if firind < firlen: return 1
    elif secind < seclen: return -1
    else: return 0

AS_prefix = "AS:i:"
def get_score(r_data, paired):
    if paired: astr = r_data[0]
    else: astr = r_data
    start_pos = astr.find(AS_prefix)+len(AS_prefix)
    fin_pos = astr.find('\t', start_pos)
    return int(astr[start_pos:fin_pos])

def get_name(line):
    return line[:line.find('\t')]

def get_name_and_score(line):
    name_pos = line.find("\t")
    score_start_pos = line.find("AS:i:", name_pos)+5
    score_fin_pos = line.find("\t", score_start_pos)
    return line[:name_pos], int(line[score_start_pos:score_fin_pos])

def asc_order(name1, name2):
    return (uporotiycmp(name1, name2) <= 0)

def main():
    # Parse arguments:

    parser = argparse.ArgumentParser()
    parser.add_argument("--pat_sam", required=True, help="Reads aligned to paternal genome")
    parser.add_argument("--mat_sam", required=True, help="Reads aligned to maternal genome")
    parser.add_argument("--o", required=True, help="Output file")
    parser.add_argument("--paired", default=0, help="Flag: If reads are paired-end")
    args = parser.parse_args()
    
    paired = int(args.paired)

    # Open output_sam; Get header:

    header = subprocess.check_output("samtools view -SH "+args.mat_sam, shell=True)
    out_stream = open(args.o, "w")
    out_stream.write(header)
    out_stream.write("@RG\tID:pat\n")

    # Open input_sam:        

    mat_count , mat_only , pat_count , pat_only = 0 , 0 , 0 , 0
    equal , mat_rand , pat_rand = 0 , 0 , 0
    bad_reads = set()

    source_m = open(args.mat_sam, 'r')
    source_p = open(args.pat_sam, 'r')

    def output_read(a_read):
        if paired:
            out_stream.write(a_read[0])
            out_stream.write(a_read[1])
        else:
            out_stream.write(a_read)
    

    def blocks_generator(fhandler):
        beg_line = fhandler.readline()
        output = [beg_line]
        name, score = get_name_and_score(beg_line)
        our_prefix = name+'\t'
        for line in fhandler:
            if line.startswith(our_prefix):
                output.append(line)
            else:
                yield output, name, score
                beg_line = line
                output = [beg_line]
                name, score = get_name_and_score(beg_line)
                our_prefix = name+'\t'
        yield output, name, score

    def correct_blocks_generator(fhandler):
        if paired:
            for block, name, score in blocks_generator(fhandler):
                if (len(block) == 2): yield block, name, score
                else: bad_reads.add(name)
        else:
            for block, name, score in blocks_generator(fhandler):
                if (len(block) == 1): yield block[0], name, score
                else: bad_reads.add(name)

    # Skip header in each file:

    m_skip = int(subprocess.check_output("samtools view -SH "+args.mat_sam+" | wc -l", shell=True).strip())
    p_skip = int(subprocess.check_output("samtools view -SH "+args.pat_sam+" | wc -l", shell=True).strip())

    for i in range(m_skip):
        source_m.readline()
    for i in range(p_skip):
        source_p.readline()
    
    # Create generator objects:
    
    pgen = correct_blocks_generator(source_p)
    mgen = correct_blocks_generator(source_m)
  
    # Merge till some EOF: 
    try:
        m_read, m_read_name, m_score = mgen.next()
        p_read, p_read_name, p_score = pgen.next()
        while 1:
            if m_read_name == p_read_name:
                if m_score > p_score:
                    output_read(m_read)
                    mat_count += 1
                elif m_score < p_score:
                    output_read(p_read)
                    pat_count += 1
                else:
                    equal += 1
                    x = random.randint(0, 1)
                    if x == 1:
                        mat_rand += 1
                        output_read(m_read)
                    else:
                        pat_rand += 1
                        output_read(p_read)
                m_read, m_read_name, m_score = mgen.next()
                p_read, p_read_name, p_score = pgen.next()
            elif asc_order(p_read_name, m_read_name):
                pat_only += 1
                output_read(p_read)
                p_read, p_read_name, p_score = pgen.next()
            else:
                mat_only += 1
                output_read(m_read)
                m_read, m_read_name, m_score = mgen.next()
    except StopIteration:
        pass

    ### Write the remain part of reads:
    for m_read, m_read_name , m_score in mgen:
        mat_only += 1
        output_read(m_read)
    for p_read, p_read_name , p_score in pgen:
        pat_only += 1
        output_read(p_read)
            

    out_stream.close(); source_m.close(); source_p.close()

    print("SUMMARY")
    print("%d reads in PAT only"%(pat_only))
    print("%d reads in MAT only"%(mat_only))
    print("%d reads where PATQ > MATQ"%(pat_count))
    print("%d reads where MATQ > PATQ"%(mat_count))
    print("%d reads where MATQ = PATQ"%(equal))
    print("%d MAT randomly selected"%(mat_rand))
    print("%d PAT randomly selected"%(pat_rand))
    print("----- %s seconds -----" % (time.time() - start_time))
    print("%d BAD read names: "%(len(bad_reads)) + " , ".join(sorted(list(bad_reads))))


if __name__ == "__main__":
    #cProfile.run("main()")
    main()
