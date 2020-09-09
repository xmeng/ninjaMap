#!/usr/bin/env python3
import sys
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

shortlist = sys.argv[1]
whole_fastq = sys.argv[2]
shortlisted_fasta = sys.argv[3]

shortlist_names = set() # Use a set instead of a list
with open(shortlist, 'r') as shortlist_file:
    shortlist_names = set([line.strip() for line in shortlist_file])

with open(shortlisted_fasta, 'w') as shortlist_fa, gzip.open(whole_fastq, 'rt') as fastq:
    for title, seq, qual in FastqGeneralIterator(fastq):
        if title in shortlist_names:
            shortlist_fa.write(f'>{title}\n{seq}\n')
