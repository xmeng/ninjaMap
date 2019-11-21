#!/usr/bin/env python3


import gzip
import os
import sys
from collections import Counter, defaultdict

import pandas as pd
import pybedtools
import pysam
from Bio import SeqIO

from ninjautils.Reads import Reads
from ninjautils.Strains import Strains

sys.path.insert(0, '/mnt')

# import argparse
# import logging
# import numpy as np
# import re
# from time import perf_counter as timer

# def read_bam_file(bamfile_name, bins):
#     total_reads = 0
#     read_objects = defaultdict()
#     perfect_alignment = defaultdict(lambda: defaultdict(int))
#     discarded_reads = set()

#     # Read the BAM file
#     bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
#     for aln in bamfile.fetch(until_eof=True):
#         read_name, mates_name, read_length, template_length = extract_read_info(aln)
#         if Reads.is_perfect_alignment(aln):
#             # read_info[read_name] = (mates_name, read_length, template_length)
#             read = Reads(read_name, mates_name, read_length, template_length)
#             read_objects[read_name] = read
#             strain_name = bins[aln.reference_name] # bins[contig_name] --> strain name
#             perfect_alignment[read_name][strain_name] += 1
#         else:
#             discarded_reads.add(read_name)
        
#     total_reads = len(discarded_reads) + len(read_objects.keys())
#     bamfile.close()
#     return(total_reads, perfect_alignment, read_objects)


def choose_primary_candidate(read, mate):
    if read.mate_has_perfect_match or mate.mate_has_perfect_match :
        common_strains_list = intersection(read.mapped_strains.keys(), mate.mapped_strains.keys())
        if len(common_strains_list) == 1:
            return common_strains_list[0].name
        else:
            return None

def calculate_coverage(bamfile_name, output_dir):
    pybedtools.set_tempdir(f'{output_dir}/tmp')
    a = pybedtools.BedTool(bamfile_name)
    df = a.genome_coverage(dz = True).to_dataframe(names=['contig','pos', 'depth'])
    pybedtools.cleanup()
    return df

def check_for_fraud(votes_file, discard_pile=1):
    fraud_status = True
    df = pd.read_csv(votes_file, index_col='Read_Name', usecols=['Read_Name', 'cSingular_Vote', 'cEscrow_Vote'])
    df['Total_Votes'] = round((df['cSingular_Vote'] + df['cEscrow_Vote']), 5)
    df = df.drop(['cSingular_Vote',  'cEscrow_Vote'], axis = 1)
    votes_df = df.groupby('Read_Name').sum()
    gt_row, gt_col = votes_df[(votes_df.Total_Votes > 1.001)].shape
    mid_row, mid_col = votes_df[(votes_df.Total_Votes > 0.001) & (votes_df.Total_Votes < 0.999)].shape
    if gt_row == 0 and mid_row == 0:
        fraud_status = False

    return fraud_status

def human_time(time):
    time = abs(time)
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    time_str = format('%02d:%02d:%02d:%02d'%(day,hour,minutes,seconds))
    return time_str

def intersection(list1, list2):
    return list(set(list1) & set(list2))

# 
# def bam_is_empty(fn):
#     if os.path.getsize(fn) > 1000000:
#         return False

#     bam = pysam.Samfile(fn, check_sq=False)
#     try:
#         bam.next()
#         return False
#     except StopIteration:
#         return True

def sort_and_index(file_name, cores=4, by='coord'):
    """ Sorts and indexes a bam file by coordinates.
    """
    if by.lower() == 'coord':
        sorted_name = file_name.replace('.bam', '') + '.sortedByCoord.bam'
        # pysam sort multithreading support doesn't work
        # pysam.sort('-@',cores,'-o',sorted_name, file_name)
        pysam.sort('-o',sorted_name, file_name)
        
    elif by.lower() == 'name':
        sorted_name = file_name.replace('.bam', '') + '.sortedByName.bam'
        pysam.sort('-n','-o',sorted_name, file_name)
    else:
        raise Exception("Bam file can only be sorted by 'coord' or 'name'.")

    pysam.index(sorted_name)
    os.remove(file_name)
    return sorted_name

def get_bam_filehandle(bamfile_name, prefix, tag):
    tmp_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    filtered_bamfile = pysam.AlignmentFile(prefix + '.ninjaMap.' + tag + '.bam', "wb", template=tmp_bamfile)
    tmp_bamfile.close()
    return filtered_bamfile

##################################
# Static methods from ninjaIndex
##################################

def get_db_metadata(fasta_list, fasta_filename):
    '''
    read the fasta files
    return bins, all_strain_obj, strains_list, concatenated fasta file
    '''
    bins = defaultdict()
    all_strain_obj = defaultdict()

    with open(fasta_filename, "w") as fasta_output:
        for input_file in fasta_list:
            # fasta_sequences = SeqIO.parse(open(input_file),'fasta')
            tmp_strain_name = os.path.basename(input_file)
            strain_name = os.path.splitext(tmp_strain_name)[0]
            strain = Strains(strain_name)
            all_strain_obj[strain_name] = strain
            # strains_list.append(strain_name)
            with open(input_file, "r") as fasta_input:
                for fasta in SeqIO.parse(fasta_input,'fasta'):
                    contig_name, sequence = fasta.id, str(fasta.seq)
                    contig_len = len(sequence)
                    strain.add_contig(contig_name, contig_len)
                    bins[contig_name] = strain_name
                    Strains.total_genome_size += contig_len

                    SeqIO.write(fasta, fasta_output , "fasta")

    return (bins, all_strain_obj, fasta_filename)

def create_bin_map(all_strain_obj, binmap_file):
    binmap = open(binmap_file, 'w')
    # Header
    # strain_name, strain_wt_1, strain_wt_2, strain_wt_3, strain_score, strain_unique_bases, contig_name, contig_length
    binmap.write('Strain_Name,Strain_Weight_Andres,Strain_Weight_Original, Strain_Weight_UScore,Strain_Uniqueness_Score,Strain_Absolute_Unique_Bases,Contig_Name,Contig_Length\n')
    for name, strain in all_strain_obj.items():
        line = ''
        first_half = strain.name +','+ \
            str(strain.adj_primary_wt) +','+ \
            str(strain.adj_primary_wt_2) +','+ \
            str(strain.adj_primary_wt_3) +','+ \
            str(strain.uniqueness_score) +','+ \
            str(strain.uniquely_covered_bases)
        for contig_name, contig_len in strain.contigs.items():
            line += str(first_half)  +','+ \
                    contig_name +','+ \
                    str(contig_len) +'\n'
        binmap.write(line)
    binmap.close()

##################################
# 01_filter_alignments.py
##################################
def filter_bam(bamfile_name, min_id, min_aln_len,
                bins, all_strain_obj, prefix):
    perfect_alignment = defaultdict(int)
    read_info = defaultdict(int)
    total_alignments = 0
    total_perfect_alignments = 0
    
    # Write intermediate files
    filtered_bamfile = get_bam_filehandle(bamfile_name, prefix, 'filtered')
    aln_stats_file = prefix +'.ninjaMap.aln_stats.csv.gz'
    aln_stats = gzip.open(aln_stats_file, mode = 'wt')
    aln_stats.write("read_name,read_length,read_mate,contig_name,contig_length,strain_name,genome_len,percent_id,perc_aln\n")
    
    # Read the BAM file
    bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    for aln in bamfile.fetch(until_eof=True):
        read_name, mates_name, read_length = Reads.extract_read_info(aln)
        read_info[read_name] += 1
        total_alignments += 1
        pId, pAln = Reads.get_aln_quality(aln)
        if Reads.is_perfect_alignment(pId, pAln, min_id, min_aln_len):
            contig_name = aln.reference_name
            strain_name = bins[contig_name] # bins[contig_name] --> strain name
            strain = all_strain_obj[strain_name]
            contig_length = strain.contigs[contig_name]
            genome_size = strain.genome_size
            # write to a new BAM file
            filtered_bamfile.write(aln)
            # write a read-strain stats data frame
                # read_name, read_length, read_mate, contig_name, contig_length, strain_name, genome_len, percent_id, perc_aln
            aln_stats.write(f'{read_name},{read_length},{mates_name},{contig_name},{contig_length},{strain_name},{genome_size},{pId},{pAln}\n')
            perfect_alignment[read_name] += 1
            total_perfect_alignments += 1

    bamfile.close()
    aln_stats.close()
    return (len(perfect_alignment.keys()), len(read_info.keys()), total_alignments, total_perfect_alignments)

##################################
# 02_filter_alignments.py
##################################
