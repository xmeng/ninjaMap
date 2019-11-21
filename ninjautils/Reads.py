#!/usr/bin/env python3

import argparse
import gzip
import logging
import numpy as np
import os
import pybedtools
import pysam
import pandas as pd
import re
import sys

from collections import defaultdict, Counter
from time import perf_counter as timer

class Reads:
    total_reads_aligned = 0
    reads_w_perfect_alignments = 0
    total_singular_reads = 0
    total_escrow_reads_kept = 0
    total_escrow_reads_discarded = 0
    total_escrow_reads = 0
    total_singular_reads_after_recruitment = 0

    def __init__(self, name, mate_name, read_length):
        self.name = name
        self.unique_name = name
        self.mates_unique_name = mate_name
        self.read_length = abs(read_length)

        self.cum_vote = 0
        self.has_voted = False
        self.in_singular_bin = False
        self.mate_has_perfect_match = False

        self.mapped_strains = defaultdict()

    def __hash__(self):
        return hash(str(self.unique_name))

    def __eq__(self, other):
        return self.unique_name == other.unique_name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def add_exact_match(self, strain):
        self.mapped_strains[strain] = strain.name
    
    def put_pair_in_singular_bin(self, mate):
        self.in_singular_bin = True
        mate.in_singular_bin = True

    def add_vote(self, vote_value):
        # vote_value = round(vote_value, 7)
        self.cum_vote += vote_value
        self.has_voted = True

    def is_fraud(self):
        '''
        for each object of this class, return True if cumulative votes > 1
        '''
        fraud = True
        # approx 1
        if (self.cum_vote < 1.001) and (self.cum_vote > 0.999):
            fraud = False

        # approx 0
        if (self.cum_vote < 0.001):
            fraud = False

        return fraud

    def get_voting_details(self, approved_strain_list):
        '''
        returns a list of 3 element lists, each containing: strain_name, singular vote value and escrow vote value
        '''
        vote_list = list()
        for strain in approved_strain_list:
            strain_name = ''
            escrow_votes = 0
            singular_votes = 0
            cumulative_vote = 0
            
            if self.cum_vote is not None:
                cumulative_vote = self.cum_vote

            if strain.name is not None:
                strain_name = strain.name

            if self.unique_name in strain.escrow_bin.keys():
                escrow_votes = strain.escrow_bin[self.unique_name]

            if self.unique_name in strain.singular_bin.keys():
                singular_votes = strain.singular_bin[self.unique_name]

            vote_list.append([strain_name, singular_votes, escrow_votes, cumulative_vote])
        return vote_list

    @staticmethod
    def get_aln_quality(aln):
        if not aln.has_tag('NM'):
            return None
        
        edit_dist = aln.get_tag('NM')
        query_len = aln.query_length
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        
        if not query_len:
            return None

        pid = (query_len - edit_dist)*100/query_len
        aln_len = aln.get_overlap(ref_start, ref_end)*100/query_len

        return (pid, aln_len)

    @staticmethod
    def is_perfect_alignment(pid, aln_len, min_perc_id = 100, min_perc_aln = 100):
        # (pid, aln_len) = get_aln_quality(aln)
        # https://www.biostars.org/p/106126/
        # return ((edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end)))
        if (pid is None):
            return None
        
        if (aln_len is None):
            return None

        return (min_perc_id <= pid <= 100) and (min_perc_aln <= aln_len <= 100)

    @staticmethod
    def parse_read_name(aln):
        '''
        Accept: AlignmentFile object from PySam
        if read name has a '/', this is the old format. 
        strip the content after the '/', return remaining
        else, return it as is.
        '''
        try:
            key, value = aln.query_name.split("/")
        except ValueError:
            return str(aln.query_name)
        else:
            return str(key)
        
    @staticmethod
    def get_unique_read_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'fwd'
        else:
            orientation =  'rev'
            
        return Reads.parse_read_name(aln) +'__'+ orientation

    @staticmethod
    def get_unique_mate_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'rev'
        else:
            orientation =  'fwd'
            
        return Reads.parse_read_name(aln) +'__'+ orientation

    @staticmethod
    def extract_read_info(aln):
        read_name = Reads.get_unique_read_name(aln)
        mate_name = Reads.get_unique_mate_name(aln)
        read_length = aln.reference_length

        return (read_name, mate_name, read_length)

class SingularShort(Reads):
    '''
    A short read pair that aligns within given parameters to a single strain
    '''
    pass

class SingularLong(Reads):
    '''
    A long read that aligns within given parameters to a single strain
    '''
    pass

class Escrow:
    '''
    A short/long read that aligns within given parameters to more than one strain
    A read (pair) CANNOT be both singular and escrow at the same time.
    '''
    pass