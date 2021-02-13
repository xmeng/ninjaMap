# input parameter: ouptut dir & input fastq dir

project=210209_MITI-MiSeq-009_JGM2J

mkdir -p /data/ampliseq/$project
cd /data/ampliseq/$project

Rscript ~/scripts/16S_R/16S.R /data/ampliseq/$project /data/myBaseSpace/$project

: <<'COMMENT'
awk 'BEGIN { FS=OFS="\t" }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}' DADA2_summary/ASVs_counts.tsv > ASVs_counts_transpose.tsv 

COMMENT

#python fix_summary.py ./DADA2_summary/Phylum_summary.txt ./DADA2_summary/Sample_stats.tsv 0
#python fix_summary.py ./DADA2_summary/Class_summary.txt ./DADA2_summary/Sample_stats.tsv 0
#python fix_summary.py ./DADA2_summary/Order_summary.txt ./DADA2_summary/Sample_stats.tsv 0
#python fix_summary.py ./DADA2_summary/Family_summary.txt ./DADA2_summary/Sample_stats.tsv 0
#python fix_summary.py ./DADA2_summary/Genus_summary.txt ./DADA2_summary/Sample_stats.tsv 10


python ~/scripts/16S_R/fix_summary.py ./Sample_stats.tsv ./ASVs_taxonomy_dada2.tsv ./ASVs_counts.tsv
