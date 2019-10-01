#!/usr/bin/env python

import pandas as pd
import re
import os
import urllib
import gzip
import tempfile
import sys

from Bio import SeqIO

## Functions
def process_raw_seedfile(seedfile):
    df = pd.read_csv(seedfile).dropna(axis=0, subset=['Strain_Name', 'FTP_link'])
    # replace all punctuations, special characters and spaces with '-'
    df['Strain_Name'] = df['Strain_Name'].apply(lambda x: re.sub('[^A-Za-z0-9]+', '-', x))
    # remove '-' if it appears at the end of the string
    df['Strain_Name'] = df['Strain_Name'].apply(lambda x: re.sub('\-+$', '', x))

    # get reference IDs from the FTP path.
    # This makes the script very specific to the current NCBI FTP paths.
    df['refIDs'] = df['FTP_link'].apply(lambda x: os.path.basename(os.path.dirname(x)))

    return df

def download_genome(genome_name, reference_id, ftp_link, temp_dir, attempt_num = 1):
    max_retries = 5
    localname = f'{temp_dir}/{genome_name}__{reference_id}.fna.gz'
    # import from ftp
    try:
        local_filename, headers = urllib.request.urlretrieve(ftp_link, localname)
    except: # catch *all* exceptions
        e = sys.exc_info()[0]
        print(f'[ERROR] Unable to retrieve reference from URL {ftp_link}\nThe following error occurred: {e.message}')
        sys.exit(1)
    else:
        outfasta, outbin, outmeta = __standardize_sequence_names(name = genome_name, 
                                                                reference_id = reference_id, 
                                                                localfile = local_filename,
                                                                temp_dir = temp_dir)

    print(f'Downloaded {genome_name}')
    return (outfasta, outbin, outmeta)

def __standardize_sequence_names(name, reference_id, localfile, temp_dir):
    output_base = f'{temp_dir}/{name}'
    
    fasta_name = f'{output_base}.fna'
    bin_name = f'{output_base}.bin.tsv'
    meta_name = f'{output_base}.metadata.tsv'
    
    outfasta = open(fasta_name,'w')
    outbin = open(bin_name,'w')
    outmeta = open(meta_name,'w')
    
    outmeta.write('AssemblyID\tGenomeName\tNewSeqHeaders\tOriginalSeqHeaders\tSeqLen\tNum_Ns\n')
    
    try: 
        with gzip.open(localfile, "rt") as handle:
            for index, record in enumerate(SeqIO.parse(handle, "fasta")):
                header = f'{name}_Node_{index}'
                # Metadata
                # Cols (tab-delim): AssemblyID, GenomeName, NewSeqHeaders, OriginalSeqHeaders, SeqLen, Num_Ns
                outmeta.write(f'{reference_id}\t{name}\t{header}\t{record.id}\t{len(record.seq)}\t{record.seq.count("N")}\n')
                outbin.write(f'{header}\t{name}\n') # contig_name\tstrain_name
                outfasta.write(f'>{header}\n{record.seq}\n') # reference fasta
    except:
        print(f'[ERROR] Could not process fasta file for {name} from {output_base}')
    else:
        os.remove(localfile)
    finally:
        outfasta.close()
        outbin.close()
        outmeta.close()

    return (fasta_name, bin_name, meta_name)

def concatenate_files(output_file, file_list, has_header):
    header_saved=False
    with open(output_file,'w') as outfile:
        for filename in file_list:
            with open(filename) as infile:
                if has_header:
                    header = next(infile)
                    if not header_saved:
                        outfile.write(header)
                        header_saved = True
                for line in infile:
                    outfile.write(line)
            os.remove(filename)

if __name__ == '__main__':
    output_dir = '/Users/sunit.jain/Research/SyntheticCommunities/NinjaMap/Database/201909/RefGenomes/20190911'
    db_name = '20190911_scv2'
    # db_ref_file example (see sheet "v2"): 
    # https://docs.google.com/spreadsheets/d/1rZUfisKKmPfLkZP1eEsNFyCKUyMm6nHQQAhNxnJJYXE/edit#gid=0
    # The file needs to have the following 2 columns: "Strain_Name" and "FTP_link"
    db_ref_file = "/Users/sunit.jain/Research/SyntheticCommunities/NinjaMap/Database/201909/RefGenomes/20190911_added_in_v2.csv"

    temp_dir = tempfile.mkdtemp()
    genomes_df = pd.DataFrame()
    db_files = pd.DataFrame()

    try:
        genomes_df = process_raw_seedfile(db_ref_file)
        db_files = genomes_df.apply(lambda row: download_genome(genome_name = row['Strain_Name'], 
                                                                reference_id = row['refIDs'], 
                                                                ftp_link = row['FTP_link'],
                                                                temp_dir = temp_dir), 
                                                                axis=1)
        db_files=db_files.apply(pd.Series).dropna()
        db_files.columns=['fastafiles', 'binfiles', 'metafiles']

        concatenate_files(output_file = f'{output_dir}/{db_name}.fna', file_list = db_files['fastafiles'], has_header=False)
        concatenate_files(output_file = f'{output_dir}/{db_name}.bin.tsv', file_list = db_files['binfiles'], has_header=False)
        concatenate_files(output_file = f'{output_dir}/{db_name}.metadata.tsv', file_list = db_files['metafiles'], has_header=True)
    except:
        e = sys.exc_info()[0]
        print(f'[FATAL] {e.message}')
        sys.exit(1)
    finally:
        os.rmdir(temp_dir)