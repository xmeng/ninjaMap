#!/usr/bin/env python3

import pysam
import argparse

p = argparse.ArgumentParser(   
    formatter_class=argparse.RawTextHelpFormatter,
    add_help=True,
    usage=argparse.SUPPRESS,
    description="""Description:
This script will allows you to filter a bam file based on certain flags.
""",
    epilog="""Examples:
python bamParser.py -bam sample.sortedByCoord.bam -id 70 -aln_len 80 -out filtered.bam
""")
# Required
p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                help='coord sorted bam file.')
p.add_argument('-id', dest='min_id', action='store', type=float, required = True,
                help='minimum percent nucleotide identity (inclusive).', default = 70)
p.add_argument('-aln_len', dest='min_aln_len', action='store', type=float, required = True,
                help='minimum percent read alignment length (inclusive).', default = 80)
p.add_argument('-se', dest='is_single_ended', action='store_true', default=False,
                help='Single ended reads were used to create this bam file; default = PE')
p.add_argument('-out', dest='output', action='store', type=str, required = True,
                help='output bam file')

args = vars(p.parse_args())
bamfile_name = args['bamfile']
min_id = args['min_id']
min_aln_len = args['min_aln_len']
out_bamfile_name = args['output']
paired = not args['is_single_ended']

def acceptable_alignment(aln, min_pid, min_paln):
    edit_dist = dict(aln.tags)['NM']
    query_len = aln.query_length
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    
    pid = (query_len - edit_dist)*100/query_len
    aln_len = aln.get_overlap(ref_start, ref_end)*100/query_len

    # https://www.biostars.org/p/106126/
    return ((pid >= min_pid) and (aln_len >= min_paln))

bam_in = pysam.AlignmentFile(bamfile_name, mode = 'rb')
bam_out = pysam.AlignmentFile(out_bamfile_name, "wb", template=bam_in)

for aln in bam_in.fetch(until_eof=True):
    if acceptable_alignment(aln, min_id, min_aln_len):
        bam_out.write(aln)

bam_out.close()
bam_in.close()