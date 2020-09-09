#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import gzip
import s3fs
import boto3

from io import BytesIO, TextIOWrapper
from collections import Counter, defaultdict

def get_genome_name(query, name_map):
    genome_name = ''
    try:
        genome_name = name_map.loc[query].Strain_Name
    except KeyError as k:
        names = query.split('_')
        genome_name = "_".join(names[:-1])
    finally:
        return genome_name

def set_bitscore_threshold(bit_list, percent = 5):
    perc = percent/100
    return max(bit_list)-max(bit_list * perc)

def generate_set_df(top_hits_df, base_col = 'UHGG', overlap_col = 'qseqid', sort_asc = True):
    df = pd.DataFrame(top_hits_df
                      .set_index(base_col)
                      .groupby(base_col)[[overlap_col]]
                      .agg(set)).reset_index()
    df['set_size'] = df[overlap_col].map(len)
    
    df = df.sort_values('set_size',ascending = sort_asc).reset_index(drop = True)
    return df

def get_confidence(top_hits_df, base_col = 'UHGG', overlap_col = 'qseqid', sort_asc = True):
    df = pd.DataFrame(top_hits_df
                      .set_index(base_col)
                      .groupby(base_col)[[overlap_col]]
                      .agg(set)).reset_index()
    df['set_size'] = df[overlap_col].map(len)
    df['conf'] = df['set_size'].apply(lambda x: 1/x)
    
    df = df.sort_values('conf',ascending = sort_asc).reset_index(drop = True)
    return df

def get_genome_confidence(df): 
    genome_conf_df = {
        'Genome' : [],
        'Conf' : [],
        'Read_Support' : [],
    }
    genome_confidence = defaultdict(float)
    genome_reads = defaultdict(int)
    
    read_confidence_df = get_confidence(df, base_col = 'qseqid', overlap_col = 'UHGG', sort_asc = False)
    
    for row_idx, row in read_confidence_df.iterrows():
        conf = row['conf']
        for genome in row['UHGG']:
            genome_confidence[genome] += conf
            genome_reads[genome] += 1

    for g,c in genome_confidence.items():
        genome_conf_df['Genome'].append(g)
        genome_conf_df['Conf'].append(c)
        genome_conf_df['Read_Support'].append(genome_reads[g])
    
    return pd.DataFrame.from_dict(genome_conf_df).sort_values('Conf',ascending = False).reset_index(drop=True)

def get_genomeset_overlaps(df, too_many = 50):
    '''
    given a list of sets
    cluster sets that are a complete subset of another
    '''
    set_dict = defaultdict(Counter)
    pid_dict = defaultdict(set)
    qcov_dict = defaultdict(set)
    registered_child = set()
    voted = set()
    read_with_too_many_matches = 0
    for parent_idx, parent_row in df.iterrows():
        parent_set = frozenset(parent_row['UHGG'])
        
        if parent_row['set_size'] >= too_many:
            parent_read_name = parent_row['qseqid']
            read_with_too_many_matches += 1
            voted.add(parent_read_name)
            continue

        # don't let something thats already a child be a parent.
        if parent_set in registered_child:
            continue

        for child_idx, child_row in df.iterrows():
            child_set = frozenset(child_row['UHGG'])
            read_name = child_row['qseqid']
            pid_set = child_row['pident']
            qcovs_set = child_row['qcovs']
            # Since the creation of a set is based entirely on the read, 
            # 1 read cannot vote for multiple sets.
            # don't allow a read to vote again
            if read_name in voted:
                continue

            if parent_set.issuperset(child_set):
                set_dict[parent_set][read_name] += 1
                pid_dict[parent_set].update(pid_set)
                qcov_dict[parent_set].update(qcovs_set)
                registered_child.add(child_set)
                voted.add(read_name)
    print(f'Reads with too many matches: {read_with_too_many_matches}')
    return set_dict, pid_dict, qcov_dict

def get_SC_name(accession_list, name_map, col = 'Genome'):
    names = set()
    for acc in accession_list:
        if acc in name_map.index:
            names.add(name_map.loc[acc][col])
    return sorted(list(names))

def get_lineage(accession_list, name_map, col = 'Lineage'):
    names = set()
    for acc in accession_list:
        if acc in name_map.index:
            names.add(name_map.loc[acc][col])
    return sorted(list(names))

def save_df(df, s3_output_dir, local_output_dir, sample_name_dir, sample_name, suffix):
    file_path = f'{sample_name_dir}/{sample_name}.{suffix}.csv'
    # to s3
    df.to_csv(f'{s3_output_dir}/{file_path}', 
                        index = False)
    # to local
    df.to_csv(f'{local_output_dir}/{file_path}',
                        index = False)

if __name__ == '__main__':
    fs = s3fs.S3FileSystem(anon=False)
    s3_blast_m8_output = sys.argv[1]
    sample_name_dir = s3_blast_m8_output.split('/')[-2]
    sample_prefix = os.path.basename(os.path.splitext(s3_blast_m8_output)[0])
    sample_name = f"{sample_name_dir}__{sample_prefix}"
    local_output_dir = 'parsed_results'
    os.makedirs(f'{local_output_dir}/{sample_name_dir}', exist_ok=True)
    s3_output_dir = 's3://czbiohub-microbiome/Sunit_Jain/scratch/NinjaMap/UHGG/NCBI_Gut_extSCv1v2/BlastN/parsed_results'
    # Metadata Files
    s3_metadata_path = 's3://czbiohub-microbiome/Sunit_Jain/scratch/NinjaMap/UHGG/NCBI_Gut_extSCv1v2/metadata_files'
    uhgg_scv1_scv2_map_file = f'{s3_metadata_path}/UHGG_SCv1_SCv2__binmap.tsv'
    genomes_metadata_all_file = f'{s3_metadata_path}/genomes_metadata.tsv'
    
    uhgg_scv1_scv2_map = pd.read_table(uhgg_scv1_scv2_map_file,
                                    usecols = ['Genome', 'UHGG', 'Genome_Size', 'Num_Contigs'],
                                    index_col = 'UHGG')
    # uhgg_scv1_scv2_map.head()

    scv1_binmap = pd.read_csv('s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/20180725/scv1/db/uniform10x_ninjaIndex.ninjaIndex.binmap.csv',
                            usecols = ['Contig_Name', 'Strain_Name'])
    scv1_binmap = scv1_binmap[['Contig_Name', 'Strain_Name']]
    # scv1_binmap.head()

    scv2_binmap = pd.read_table('s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Narrow/20190911/scv2/db/20190911_scv2.bin.tsv',
                            header = None,
                            names = ['Contig_Name', 'Strain_Name'])
    # scv2_binmap.head()
    genomes_metadata_all = pd.read_table(genomes_metadata_all_file,
                                        index_col = 'Genome',
                                        usecols = ['Genome','Genome_type','Lineage'],
                                    low_memory=False)

    genomes_metadata = genomes_metadata_all[genomes_metadata_all.Genome_type == 'Isolate']
    del genomes_metadata_all
    
    contig_strain_map = pd.concat([scv1_binmap, scv2_binmap]).drop_duplicates()
    contig_strain_map['Contig_Name'].value_counts(ascending = False)
    contig_strain_map.set_index('Contig_Name', inplace = True)

    blast_hits = pd.read_table(s3_blast_m8_output,
                            names = ['qseqid','sseqid','pident','length',
                                    'mismatch','gapopen','qstart','qend',
                                    'sstart','send','evalue','bitscore',
                                    'qlen','slen','qcovs'],
                            index_col = 'qseqid')
    blast_hits['UHGG'] = blast_hits['sseqid'].apply(get_genome_name, name_map=contig_strain_map)
    blast_hits['bit_threshold'] = blast_hits.groupby('qseqid')['bitscore'].apply(set_bitscore_threshold)
    blast_hits['best_bitscore'] = blast_hits.groupby('qseqid')['bitscore'].apply(set_bitscore_threshold, percent = 0)
    # print(blast_hits.sample(5))

    blast_hits_meta = blast_hits.join(other = uhgg_scv1_scv2_map, on='UHGG', how = 'left')

    blast_hits_wLineage = pd.merge(blast_hits_meta.reset_index(), genomes_metadata, left_on='UHGG', right_on='Genome', how = 'left')
    save_df(blast_hits_wLineage,s3_output_dir,local_output_dir, sample_name_dir, sample_name, suffix='wLineage')
    # blast_hits_wLineage.head()

    blast_hits_top = blast_hits_wLineage.query('bitscore >= best_bitscore')
    save_df(blast_hits_top,s3_output_dir,local_output_dir, sample_name_dir, sample_name, suffix='wLineage.bitScore_top')
    # blast_hits_top.set_index('qseqid', inplace = True)
    # blast_hits_top.head()
    genome_confidence_df = get_genome_confidence(blast_hits_top)
    genome_confidence_meta_tmp = pd.merge(genome_confidence_df.reset_index(drop=True), genomes_metadata, left_on='Genome', right_on='Genome', how = 'left')
    genome_confidence_meta = pd.merge(genome_confidence_meta_tmp.reset_index(drop=True), uhgg_scv1_scv2_map, left_on='Genome', right_on='UHGG', how = 'left')
    save_df(genome_confidence_meta,s3_output_dir,local_output_dir,sample_name_dir, sample_name, suffix='wLineage.genome_confidence')

    genome_set_df = generate_set_df(blast_hits_top, base_col = 'qseqid', overlap_col = 'UHGG', sort_asc=False)
    pident_set_df = generate_set_df(blast_hits_top, base_col = 'qseqid', overlap_col = 'pident', sort_asc=False).drop(['set_size'], axis=1)
    qcovs_set_df = generate_set_df(blast_hits_top, base_col = 'qseqid', overlap_col = 'qcovs', sort_asc=False).drop(['set_size'], axis=1)

    master_df = genome_set_df.merge(pident_set_df,on='qseqid').merge(qcovs_set_df,on='qseqid')

    genome_set_dict, genome_pid_dict, genome_qcov_dict = get_genomeset_overlaps(master_df)
    
    output_df = {
        'Possible_Strains' : [],
        'Cluster_Members' : [],
        'Num_Reads' : [],
        'PiD' : [],
        'QCovs' : [],
        'Reads' : [],
        'SC_Name' : [],
        'UHGG_Lineage' : []
    }
    for key1, key_dict in genome_set_dict.items():
        query = "|".join(key1)
        num_members = len(key1)
        query_SC_name = "|".join(get_SC_name(key1,uhgg_scv1_scv2_map))
        query_lineage = "|".join(get_lineage(key1,genomes_metadata))
        pids = "|".join([str(x) for x in genome_pid_dict[key1]])
        qcovs = "|".join([str(x) for x in genome_qcov_dict[key1]])
        num_reads = 0
        read_set = set()
        for key2, value in genome_set_dict[key1].items():
            num_reads += value
            read_set.add(str(key2))
        read_set_text = "|".join(read_set)

        output_df['Possible_Strains'].append(query)
        output_df['Cluster_Members'].append(num_members)
        output_df['Num_Reads'].append(num_reads)
        output_df['PiD'].append(pids)
        output_df['QCovs'].append(qcovs)
        output_df['Reads'].append(read_set_text)
        output_df['SC_Name'].append(query_SC_name)
        output_df['UHGG_Lineage'].append(query_lineage)
            # fh.write(f'{query}\t{num_members}\t{num_reads}\t{pids}\t{qcovs}\t{read_set_text}\t{query_SC_name}\t{query_lineage}\n')

    save_df(pd.DataFrame.from_dict(output_df),s3_output_dir,local_output_dir,sample_name_dir, sample_name, suffix='wLineage.bitScore_top.genome_sets')