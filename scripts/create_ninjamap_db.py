#!/usr/bin/env python

import pandas as pd
import re
import os
import urllib
import gzip

from Bio import SeqIO

output_dir='.'
# genome_links_file = "/Users/sunit.jain/Research/SyntheticCommunities/NinjaMap/Database/201909/RefGenomes/20190911_added_in_v2.csv"

data = {"provided_names" : ["Adlercreutzia equolifaciens DSM 19450",  "Alistipes finegoldii DSM 17242",
    "Alistipes indistinctus YIT 12060/DSM 22520",  "Alistipes onderdonkii DSM 19147",
    "Alistipes shahii WAL 8301/DSM 19121",  "Bacteroides rodentium DSM 26882",
    "Clostridiales bacterium VE202-03",  "Clostridiales bacterium VE202-14",
    "Oscillibacter sp. KLE 1728",  "Subdoligranulum sp. 4_3_54A2FAA",
    "Alistipes senegalensis JC50/DSM 25460",  "Alistipes ihumii AP11",
    "Bilophila wadsworthia ATCC 49260",  "Blautia wexlerae DSM 19850",
    "Burkholderiales bacterium 1_1_47",  "Butyricimonas virosa DSM 23226",
    "Clostridiales bacterium VE202-27",  "Clostridium sp. ATCC 29733 VPI C48-50",
    "Intestinimonas butyriciproducens DSM 26588",  "Odoribacter splanchnicus DSM 20712",
    "Ruminococcus gauvreauii DSM 19829",  "Blautia sp. KLE 1732 (HM 1032)"],
  "ftp_links" : ["ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/478/885/GCF_000478885.1_ASM47888v1/GCF_000478885.1_ASM47888v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/265/365/GCF_000265365.1_ASM26536v1/GCF_000265365.1_ASM26536v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/231/275/GCF_000231275.1_Alis_indi_YIT_V1/GCF_000231275.1_Alis_indi_YIT_V1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/374/505/GCF_000374505.1_ASM37450v1/GCF_000374505.1_ASM37450v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/575/GCF_000210575.1_ASM21057v1/GCF_000210575.1_ASM21057v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/614/125/GCF_000614125.1_ASM61412v1/GCF_000614125.1_ASM61412v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/508/885/GCF_000508885.1_ASM50888v1/GCF_000508885.1_ASM50888v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/509/025/GCF_000509025.1_ASM50902v1/GCF_000509025.1_ASM50902v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/469/425/GCF_000469425.1_ASM46942v1/GCF_000469425.1_ASM46942v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/635/GCF_000238635.1_Subdol_sp_4_3_54A2FAA_V1/GCF_000238635.1_Subdol_sp_4_3_54A2FAA_V1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/312/145/GCF_000312145.1_ASM31214v1/GCF_000312145.1_ASM31214v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/321/205/GCF_000321205.1_ASM32120v1/GCF_000321205.1_ASM32120v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/701/705/GCF_000701705.1_ASM70170v1/GCF_000701705.1_ASM70170v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/484/655/GCF_000484655.1_ASM48465v1/GCF_000484655.1_ASM48465v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/144/975/GCF_000144975.1_ASM14497v1/GCF_000144975.1_ASM14497v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/519/105/GCF_000519105.1_ASM51910v1/GCF_000519105.1_ASM51910v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/509/145/GCF_000509145.1_ASM50914v1/GCF_000509145.1_ASM50914v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/466/605/GCF_000466605.1_ASM46660v1/GCF_000466605.1_ASM46660v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/096/335/GCF_003096335.1_ASM309633v1/GCF_003096335.1_ASM309633v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/190/535/GCF_000190535.1_ASM19053v1/GCF_000190535.1_ASM19053v1_genomic.fna.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/425/525/GCF_000425525.1_ASM42552v1/GCF_000425525.1_ASM42552v1_genomic.fna.gz",    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/466/565/GCF_000466565.1_ASM46656v1/GCF_000466565.1_ASM46656v1_genomic.fna.gz"
      ]}

# genomes_df = pd.DataFrame(data)
# genomes_df.head()

## Functions
def process_raw_seedfile(seedfile):
    df = pd.DataFrame(seedfile)
    # replace all punctuations, special characters and spaces with '-'
    df['provided_names'] = df['provided_names'].apply(lambda x: re.sub('[^A-Za-z0-9]+', '-', x))
    # remove '-' if it appears at the end of the string
    df['provided_names'] = df['provided_names'].apply(lambda x: re.sub('\-+$', '', x))

    # get reference IDs from the FTP path.
    # This makes the script very specific to the current NCBI FTP paths.
    df['refIDs'] = df['ftp_links'].apply(lambda x: os.path.basename(os.path.dirname(x)))

    return df

def download_genome(genome_name, reference_id, ftp_link):
    localname = f'{output_dir}/{genome_name}__{reference_id}.fna.gz'
    # import from ftp
    try:
        local_filename, headers = urllib.request.urlretrieve(ftp_link, localname)
    except:
        print(f'[ERROR] Unable to retrieve reference from URL: {ftp_link}')
        return ''
    else:
        __standardize_sequence_names(name = genome_name, reference_id = reference_id, localfile = local_filename)

    print(f'Downloaded {genome_name}')
    return local_filename

def __standardize_sequence_names(name, reference_id, localfile):
    output_base = f'{output_dir}/{name}'
    outfasta = open(f'{output_base}.fna','w')
    outbin = open(f'{output_base}.bin.tsv','w')
    outmeta = open(f'{output_base}.metadata.tsv','w')
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

if __name__ == '__main__':
    genomes_df = process_raw_seedfile(data)
    # genomes_df.shape
    genomes_df['local_names'] = genomes_df.apply(lambda row: download_genome(genome_name = row['provided_names'], 
                                                                            reference_id = row['refIDs'], 
                                                                            ftp_link = row['ftp_links']), 
                                                            axis=1)
    # genomes_df.apply(lambda row: standardize_sequence_names(name = row["provided_names"], 
    #                                                         localfile = row["local_names"]),
    #                             axis=1)
