#!/usr/bin/env python3

import gzip
import os
import sys
import pysam
import subprocess

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

def sort_and_index2(file_name, by='coord', cores=4, memPerCore=6):
    sorted_name = None
    if by.lower() == 'coord':
        sorted_name = file_name.replace('.bam', '') + '.sortedByCoord.bam'

        sort_command = f'samtools sort -@ {cores} -m {memPerCore}G -o {sorted_name} {file_name}'
        subprocess.run(sort_command, shell=True, check=True)
        
        subprocess.run(f'samtools index -@ {cores} {sorted_name}', shell=True, check=True)
    elif by.lower() == 'name':
        sorted_name = file_name.replace('.bam', '') + '.sortedByName.bam'
        
        sort_command = f'samtools sort -n -@ {cores} -m {memPerCore}G -o {sorted_name} {file_name}'
        subprocess.run(sort_command, shell=True, check=True)
    else:
        raise Exception("Bam file can only be sorted by 'coord' or 'name'.")

    return sorted_name

def bam2sam(bamfile_name):
    samfile_name = bamfile_name.replace('.bam', '.sam')
    subprocess.run(f'samtools view -h -@ {cores} -o {bamfile_name} {samfile_name}', shell=True, check=True)
    return samfile_name

def sam2bam(samfile_name, sort_by='coord'):
    bamfile_name = samfile_name.replace('.sam', '.bam')
    subprocess.run(f'samtools view -h -b -@ {cores} -o {samfile_name} {bamfile_name}', shell=True, check=True)
    sorted_bam = sort_and_index(bamfile_name, by=sort_by)
    return sorted_bam

def split_bamfile(bamfile_name, primary_df, escrow_df, 
                prefix, by='coord', cores=4, memPerCore=6):
    samfile = bam2sam(bamfile_name)
    
    primary_sam = f'{prefix}.ninjaMap.primary.sam'
    primary_samhandle = open(primary_sam, 'w')

    escrow_sam = f'{prefix}.ninjaMap.escrow.sam'
    escrow_samhandle = open(escrow_sam, 'w')
    
    with open(samfile, "r") as samhandle:
        samline = samhandle.readline()
        columns = samline.split('\t')

        while line:
            if len(columns) < 5:
                primary_samhandle.write(line)
                escrow_samhandle.write(line)
            elif columns[0] in primary_df.read_names.to_numpy():
                primary_samhandle.write(line)
            elif columns[0] in escrow_df.read_names.to_numpy():
                escrow_samhandle.write(line)

    primary_bam = sam2bam(primary_sam, sort_by=by)
    escrow_bam = sam2bam(escrow_sam, sort_by=by)

    return (primary_bam, escrow_bam)

    # read sam file line by line:
    #     if line has less than 5 columns:
    #         write to both files
    #     elif columns > 5:
    #         if col[0] in primary_df.read_names:
    #             write to primary sam
    #         elif col[0] in escrow_df.read_names:
    #             write to escrow sam
    # sam2bam(primary.sam, sort_by='coord')
    # sam2bam(escrow.sam, sort_by='coord')

sorted_name = sort_and_index2(sys.argv[1], by = 'name')
print(f'{sorted_name} was successfully sorted and indexed')
