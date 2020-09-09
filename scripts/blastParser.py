#!/usr/bin/env python3
import argparse
import gzip
import logging
import os
import pprint
import sys
from time import *

import numpy as np
import pandas as pd
import s3fs
from tqdm import tqdm


def usage():
    p = argparse.ArgumentParser(   
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
        usage=argparse.SUPPRESS,
        description="""Description:
    Usage: blastParser.py -blast blastn.m6.txt.gz -out blastn.parsed.tsv.gz -min_cov 99 -min_id 90 -bit_dev 5 -evalue 1e-6
    """)

    # Required
    p.add_argument('-blast', dest='blastfile', action='store', type=str, required = True,
                    help='Tab-delim blast output with qlen, slen qcovs and sscinames as extra columns')
    p.add_argument('-prefix', dest='prefix', action='store', type=str, required = True,
                    help='tab-delimited file with taxa assigned to each contig based on blast hit')
    p.add_argument('-bit_dev', dest='bit_dev', action='store', type=float, default = 5,
                    help='bitscore deviation allowed from top hit')
    p.add_argument('-evalue', dest='evalue', action='store', type=float, default = 1e-6,
                    help='maximum evalue')
    p.add_argument('-min_id', dest='min_id', action='store', type=float, default = 90,
                    help='minimum query percent identity')
    p.add_argument('-min_cov', dest='min_cov', action='store', type=float, default = 99,
                    help='minimum query coverage')

    return vars(p.parse_args())

def get_contig_names(query):
    names = query.split('_')
    return "_".join(names[0:-1]), names[-1]

def set_bitscore_threshold(bit_list, percent):
    perc = percent/100
    return max(bit_list)-max(bit_list * perc)


def parse_names(names_dmp_path):
    names = {}
    print(strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid names')
    wc_output = subprocess.check_output(['wc', '-l', names_dmp_path])
    wc_list = wc_output.split()
    number_of_lines = int(wc_list[0])
    with open(names_dmp_path) as names_dmp:
        for line in tqdm(names_dmp, total=number_of_lines):
            taxid, name, _, classification = line.strip('\t|\n').split('\t|\t')[:4]
            taxid = int(taxid)
            # Only add scientific name entries
            scientific = classification == 'scientific name'
            if scientific:
                # line_list[1] = line_list[1].replace(' ', '_')
                names.update({taxid:name})
    return(names)

def parse_nodes(nodes_dmp_path):
    print(strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid nodes')
    wc_output = subprocess.check_output(['wc', '-l', nodes_dmp_path])
    wc_list = wc_output.split()
    number_of_lines = int(wc_list[0])
    nodes_dmp = open(nodes_dmp_path)
    root_line = nodes_dmp.readline()
    nodes = {}
    nodes.update({1:{'parent':1, 'rank':'root'}})
    for line in tqdm(nodes_dmp, total=number_of_lines):
        child, parent, rank = line.split('\t|\t')[:3]
        parent, child = map(int,[parent, child])
        nodes.update({child:{'parent':parent, 'rank':rank}})
    nodes_dmp.close()
    return(nodes)

def isCommonAncestor(parent_taxid, child_taxid, nodes_dict):
    ancestor_taxid = child_taxid
    while ancestor_taxid != 1:
        if parent_taxid == ancestor_taxid:
            return True
        ancestor_taxid = nodes_dict[ancestor_taxid]['parent']
    return False

def get_lineage(child_taxid, nodes_dict):
    ancestor_taxid = child_taxid
    while ancestor_taxid != 1:
        ancestor_taxid = nodes_dict[ancestor_taxid]['parent']
    return False

if __name__ == '__main__':
    args = usage()

    blast_input = args['blastfile']
    output_file = f'{args["prefix"]}.csv.gz'
    best_output_file = f'{args["prefix"]}.best_guess.csv.gz'

    min_cov = args['min_cov']
    min_id = args['min_id']
    max_eval = args['evalue']
    bit_dev = args['bit_dev']

    rank_priority = [
    'species',
    'genus',
    'family',
    'order',
    'class',
    'phylum',
    'superkingdom',
    'root']

    # Read blast tab results
    blast_hits = pd.read_csv(blast_input,
                            names = ['qseqid','sseqid','pident','length',
                                    'mismatch','gapopen','qstart','qend',
                                    'sstart','send','evalue','bitscore',
                                    'qlen','slen','qcovs','sscinames'],
                            # index_col = 'qseqid',
                            sep = '\t')
    
    # Parse Contig name and gene num
    blast_hits['contig'], blast_hits['gene_num'] = zip(*blast_hits['qseqid'].map(get_contig_names))
    # Get genes per contigs
    genes_per_contig = blast_hits[['contig','gene_num']].groupby('contig')['gene_num'].nunique().reset_index()

    # Process Blast hits
    blast_hits.set_index('qseqid', inplace = True)
    blast_hits['bit_threshold'] = blast_hits.groupby(['qseqid'])['bitscore'].apply(set_bitscore_threshold, percent = bit_dev)
    blast_hits.query(f'bitscore >= bit_threshold and qcovs >= {min_cov} and pident >= {min_id} and evalue <= {max_eval}', 
                    inplace = True)
    blast_hits.reset_index(inplace=True)

    # Aggregate blast hits
    agg_df = blast_hits[['contig', 'sscinames','gene_num']].groupby(['contig', 'sscinames']).agg({'gene_num' : 'nunique'}).reset_index()

    # Calculate the Taxa assignment for each contig
    contig_conf = agg_df.join(genes_per_contig.set_index('contig'), how='left', on = 'contig', lsuffix='_assigned', rsuffix='_total')
    contig_conf['conf'] = contig_conf['gene_num_assigned']/contig_conf['gene_num_total']
    
    contig_conf.set_index('contig').sort_values(by = ['contig','conf'], ascending = [True, False]).to_csv(output_file)
    contig_conf[contig_conf.groupby(['contig'])['conf'].transform(max) == contig_conf['conf']].to_csv(best_output_file, index=False)
