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

class Parser:
    @staticmethod
    def is_perfect_alignment(aln):
        edit_dist = dict(aln.tags)['NM']
        query_len = aln.query_length
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        
        # https://www.biostars.org/p/106126/
        return ((edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end)))
    
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
            
        return Parser.parse_read_name(aln) +'__'+ orientation
    
    @staticmethod
    def get_unique_mate_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'rev'
        else:
            orientation =  'fwd'
            
        return Parser.parse_read_name(aln) +'__'+ orientation
    
    @staticmethod
    def extract_read_info(aln):
        read_name = Parser.get_unique_read_name(aln)
        mate_name = Parser.get_unique_mate_name(aln)
        read_length = aln.reference_length
        template_length = aln.template_length

        return (read_name, mate_name, read_length, template_length)
