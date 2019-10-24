#!/usr/bin/env python3

import pysam
import argparse
import sys
from collections import defaultdict

def usage():
    p = argparse.ArgumentParser(   
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
        usage=argparse.SUPPRESS,
        description="""Description:
    This script will allow you to filter a bam file based on certain flags.
    """,
        epilog="""Examples:
    python bamParser.py -bam sample.sortedByCoord.bam -id 70 -aln_len 80 -out filtered.bam
    or
    python bamParser.py -bam sample.sortedByCoord.bam -id 70 -aln_len 80 -reads readnames.txt -out filtered.bam
    or
    python bamParser.py -bam sample.sortedByCoord.bam -id 70 -aln_len 80 -reads readnames.txt -out filtered.tsv -out_fmt tsv
    """)
    # Required
    p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                    help='coord sorted bam file.')
    p.add_argument('-out', dest='output', action='store', type=str, required = True,
                    help='output file name')

    p.add_argument('-id', dest='min_id', action='store', type=float,
                    help='minimum percent nucleotide identity (inclusive).', default = 0)
    p.add_argument('-aln_len', dest='min_aln_len', action='store', type=float,
                    help='minimum percent read alignment length (inclusive).', default = 0)
    p.add_argument('-reads', dest='read_list', action='store', type=str,
                    help='only print reads mentioned in this file. 3 tab sep columns w/o header: read_name, organism_name, read_length')
    p.add_argument('-out_fmt', dest='output_fmt', action='store', type=str, default = 'bam',
                    help='output format. Choose either "tsv" or "bam". Default is bam')

    return vars(p.parse_args())

def acceptable_alignment(aln, min_pid, min_paln):
    edit_dist = dict(aln.tags)['NM']
    query_len = aln.query_length
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    
    pid = (query_len - edit_dist)*100/query_len
    aln_len = aln.get_overlap(ref_start, ref_end)*100/query_len

    # https://www.biostars.org/p/106126/
    return ((pid >= min_pid) and (aln_len >= min_paln))

def filter_bam(bamfile_name, min_id, min_aln_len,
              output_file_name, output_fmt, read_file, read_filter):
    # open the BAM file for reading
    bam_in = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    # figure out the output file handle
    if output_fmt == 'tsv':
        out_file = open(output_file_name, 'w')
    else:
        out_file = pysam.AlignmentFile(output_file_name, "wb", template=bam_in)

    # create the read filter dict
    reads = defaultdict(int)
    if read_filter:
        with open(read_file, 'r') as readFile:
            for line in readFile:
                read, subject, read_bp = line.strip().split('\t')
                reads[read] += 1

    # Go through the bam file, line-by-line till the end of file
    for aln in bam_in.fetch(until_eof=True):
        # Continue until you find a read that's in the read filter
        if read_filter and not reads[aln.query_name]:
            continue

        # Write to appropriate output 
        # if alignment matches all the minimum alignment criteria?
        if acceptable_alignment(aln, min_id, min_aln_len):
            if read_filter:
                out_file.write(f'{aln.query_name}\t{aln.reference_name}\n')
            else:
                out_file.write(aln)
    
    # close the files
    out_file.close()
    bam_in.close()

if __name__ == "__main__":
    args = usage()

    bamfile_name = args['bamfile']
    output_file_name = args['output']
    output_fmt = args['output_fmt']

    min_id = args['min_id']
    min_aln_len = args['min_aln_len']
    read_filter = False
    read_file = ''
    if args['read_list']:
        read_filter = True
        read_file = args['read_list']
    else:
        output_fmt = 'bam'

    filter_bam(bamfile_name,min_id, min_aln_len,output_file_name,output_fmt,read_file,read_filter)
    sys.exit(0)