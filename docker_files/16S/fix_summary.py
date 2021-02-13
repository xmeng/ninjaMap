import numpy as np
import pandas as pd
import sys
import os.path
from itertools import compress

# 3 input files and 2 output files
# output1: ./DADA2_summary/ASVs_taxonomy_count.tsv
# output2: ./DADA2_summary/Genus_summary.tsv

#sumfile = sys.argv[1] # input genus summary table
statfile = sys.argv[1] # input sample reads stat table
#min_count = int(sys.argv[3]) # minimum frequency count
#asv_table = sys.argv[4] # ASV table transpose
taxonomy_table = sys.argv[2] # ASVs_taxonomy_dada2.tsv
asv_table_origin = sys.argv[3] # ASV table

#basename = os.path.basename(sumfile)
#base = os.path.splitext(basename)[0]
base = "Genus_summary"
# convert feature names and sample names
#df = pd.read_csv(sumfile, sep='\t', index_col=0)
# Samples are the row names written by `DADA2`. They are of the form
# sample_name.fastq so we remove the .fastq and zero pad them because
# they are just numeric in our original mapping file.
#df.columns.values[0] = 'Sample'
#df.index = ['%s' % i.strip().split('.fastq')[0] for i in df.index]

# Get the sample stats
stat = pd.read_csv(statfile, sep='\t', index_col=0)
sel_stat_df=stat[['input', 'nonchim']]  #.set_index('sample')
total_filtered=stat[['nonchim']]
# read ASV table
#asv = pd.read_csv(asv_table, sep='\t', index_col=0)
# read taxonomy table
taxonomy = pd.read_csv(taxonomy_table, sep='\t', index_col=0)
counts = pd.read_csv(asv_table_origin, sep='\t', index_col=0)

#######################################################################
if not os.path.exists('DADA2_summary'):
        os.mkdir('DADA2_summary')
# concat  ASV taxonomy and count table
# replace empty cell with "NA"
taxonomy = taxonomy.replace(np.nan, 'NA')

tc = pd.concat([taxonomy, counts], axis=1).reset_index()

def prevalence(row):
    c=0
    for i, item in enumerate(row):
        if row[i] > 0:
            c = c + 1
    return c
ct = pd.DataFrame()
#ct['prevalence'] =  counts.apply(lambda row : prevalence(row), axis=1)
ct['prevalence'] = counts.apply(lambda row : np.count_nonzero(row), axis = 1)
#ct.to_csv("ct.txt", sep='\t')
#counts.apply(lambda x: x.count() if x > 0, axis=1)

tc_sample_list = list(counts.columns)
tc_sample_list  = ["index","ASV" , "Classification", "Merged_ASVs"] + tc_sample_list
tc.to_csv("tc.txt", sep='\t')

grouped_multiple = tc.groupby(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'], as_index=False).transform('sum')
#taxonomy['taxonomy'] = taxonomy[taxonomy.columns[0:]].apply( lambda x: ';'.join(x.astype(str)), axis=1)
taxonomy1 = taxonomy[taxonomy.columns[0:]].apply( lambda x: ';'.join(x.astype(str)), axis=1).reset_index()
ct['taxonomy'] = taxonomy[taxonomy.columns[0:]].apply( lambda x: ';'.join(x.astype(str)), axis=1) #.reset_index()
grouped_multiple = pd.concat([taxonomy1, grouped_multiple], axis=1, ignore_index=True, sort=False ).reset_index()
#taxonomy1.to_csv("./taxonomy1.txt", sep='\t')
#counts.to_csv("./counts.txt", sep='\t')
ASVs_taxonomy_count = pd.concat([ct['taxonomy'], ct['prevalence'], counts], axis=1, ignore_index=False, sort=False ).sort_values('prevalence', ascending = False).reset_index()
ASVs_taxonomy_count.to_csv("./DADA2_summary/ASVs_taxonomy_count.tsv", sep='\t')

# dropping ALL duplicte ASVs
grouped_multiple.drop_duplicates(subset=grouped_multiple.columns[2], keep='first', inplace=True)

#print(grouped_multiple)
# assingn sample name, then transpose the dataframe
grouped_multiple.columns = tc_sample_list
df=grouped_multiple[grouped_multiple.columns[2:]].set_index('Merged_ASVs').transpose()

# drop the 2nd line "Classification", make it as column header
classification_header = df.iloc[0] #grab the first row for the header.
df = df[1:] #take the data less the header row.
merged_ASVs_list = list(df.columns)
df.columns = classification_header #set the header row as the df header.
df.to_csv("grouped_multiple.txt", sep='\t')

df= df.fillna(0)
#grouped_multiple = grouped_multiple.groupby(['index'])
#grouped_multiple.to_csv("grouped_index.txt", sep='\t')
#######################################################################
# Filter low frequency counts
#df = df.iloc[1:, 0:][ df > min_count ].fillna(0).astype(int)
# Do stats
df2 = pd.DataFrame()
df2['Reads_counts_genus_level/Total_filtered_reads(%)'] = round(df.iloc[:, 0:].sum(axis=1)/total_filtered.iloc[:,0]*100,2)

#df2['missed_counts'] = df.apply(lambda row: find_miss_counts(row, total_filtered))  #, axis=1)
#df2['asv_counts/Filtered_reads_counts(%)'] = round(asv.iloc[:, 0:].sum(axis=1)/total_filtered.iloc[:,0]*100,2)
df2['1st_Read'] = df.iloc[:, 0:].max(axis=1)
df2['1st_Abundance(%)'] = round(df.iloc[:, 0:].max(axis=1)/df.iloc[:, 0:].sum(axis=1)*100,2)
# index of maximum item
df2['1st_Taxonomy'] = df.iloc[:, 0:].idxmax(axis=1)
df2['2nd_Read'] = df.iloc[:, 0:].apply(lambda row: row.nlargest(2).values[-1],axis=1)
df2['2nd_Abundance(%)'] = round(df.iloc[:, 0:].apply(lambda row: row.nlargest(2).values[-1],axis=1)/df.iloc[:, 0:].sum(axis=1)*100,2)
df2['2nd_Taxonomy'] = df.iloc[:, 0:].apply(lambda row: row.nlargest(2).idxmin(), axis=1)
df2['3rd_Read'] = df.iloc[:, 0:].apply(lambda row: row.nlargest(3).values[-1],axis=1)
df2['3rd_Abundance(%)'] = round(df.iloc[1:, 0:].apply(lambda row: row.nlargest(3).values[-1],axis=1)/df.iloc[:, 0:].sum(axis=1)*100,2)
df2['3rd_Taxonomy'] = df.iloc[:, 0:].apply(lambda row: row.nlargest(3).idxmin(), axis=1)

df2['PASS/FAIL'] = df2.iloc[:, 0:].apply(lambda x: 'PASS' if x['1st_Abundance(%)'] >=99 and x['2nd_Abundance(%)'] <= 0.5  else 'FAIL', axis=1)

# change the order of columns
cols = ['PASS/FAIL']  + [col for col in df2 if col != 'PASS/FAIL']
df2 = df2[cols]


"""

sample_dict = {}
sample_list = df2.index[df2['Reads_counts_genus_level/Total_filtered_reads(%)'] < 100 ].tolist()

#taxonomy = taxonomy.replace(np.nan, 'NA')

for sample in sample_list:
    asv_list = asv.loc[sample].apply(lambda x: x > 0).tolist()
    asv_index = list(compress(range(len(asv_list)), asv_list))
    asv_value = asv.loc[sample].tolist()
    #asv_names = asv.columns.tolist()
    out_list = []
    for i in asv_index:
        #print(asv.columns[i])
        t_list = taxonomy.loc[asv.columns[i]].tolist()
        #print(t_list)
        genus = ';'.join(map(str, t_list))

        if  genus.endswith('NA'):
            genus = genus + '(' + str(asv_value[i]) +')'
            out_list.append(genus)
            #print (sample, '\t' ,genus)
    sample_dict[sample] = ''.join(map(str, out_list))
    #print (sample_dict)
    #print(sample, '\t', ''.join(map(str, out_list)))

#convert the dic to dataframe
df_missed_counts = pd.DataFrame.from_dict(sample_dict, orient='index', columns=['Missed_Classification'])

sum_df = pd.concat([sel_stat_df, df_missed_counts, df2, df], axis=1) #, sort=False).reset_index().rename(columns={'index':'sample'})
"""

sum_df = pd.concat([sel_stat_df, df2, df], axis=1)
# Added the ASVC as the 1st row in output table
sum_df.loc[' '] = ["","", "", "", "", "", "", "", "", "", "", "", ""] + merged_ASVs_list
newIndex=[' ']+[ind for ind in sum_df.index if ind!=' ']
sum_df=sum_df.reindex(index=newIndex)
#Rename col names
sum_df.rename(columns = {"input" : "Input_reads"}, inplace = True)
sum_df.rename(columns = {"nonchim" : "Filtered_reads"}, inplace = True)

#debug
#sum_df.to_csv("test.out", sep='\t')

# reorder the columns
#cols = list(sum_df)
#cols = cols[:2] + cols[-2:] + cols[2:-2]
#sum_df = sum_df[cols]

out_name= './DADA2_summary/'+ base + '.tsv'

sum_df.to_csv(out_name, sep='\t')
