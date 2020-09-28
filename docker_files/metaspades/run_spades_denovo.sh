#!/bin/bash -x

set -euoE pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
DEFAULT=$coreNum
coreNum=${coreNum:-$DEFAULT}

if [ -z "${coreNum:-}" ]; then
  echo "coreNum was not set"
fi

# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#S3OUTPUTPATH = "${3}"

# Setup directory structure
LOCAL_DB_PATH=${LOCAL}/databases
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
#QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/SPAdes"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT_1="${LOCAL_OUTPUT}/fastqc_raw"
FASTQC_OUTPUT_2="${LOCAL_OUTPUT}/fastqc_filtered"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
QC_FASTQ="${LOCAL_OUTPUT}/trimmed_fastq"
#FASTQ_NAME=${fastq1%/*}

TMP_BWT_OUTPUT="${OUTPUTDIR}/bowtie2"
BWT_OUTPUT="${LOCAL_OUTPUT}/bowtie2"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT_1}" "${FASTQC_OUTPUT_2}" "${BAMQC_OUTPUT}"
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}" "${TMP_BWT_OUTPUT}"

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

#cp  ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#cp  ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT_1} \
"${RAW_FASTQ}/read1.fastq.gz" \
"${RAW_FASTQ}/read2.fastq.gz"

# Constant definitions for bbduk
adapterFile="adapters,phix"
trimQuality=${trimQuality:-25}
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim short reads, -eoom exits when out of memory
# Filter out reads that have an average entropy of under 0.5. A homopolymer such as GGGGGGGGGG would have entropy of zero
timem bbduk.sh -Xmx32g tbo -eoom hdist=1 qtrim=rl ktrim=r \
    entropy=0.5 entropywindow=50 entropyk=5 \
    in1="${RAW_FASTQ}/read1.fastq.gz" \
    in2="${RAW_FASTQ}/read2.fastq.gz" \
    out1="${QC_FASTQ}/read1_trimmed.fastq.gz" \
    out2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
    ref=${adapterFile} \
    k="${kmer_value}" \
    mink="${min_kmer_value}" \
    trimq="${trimQuality}" \
    minlen="${minLength}" \
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" |\
    tee -a ${LOG_DIR}/bbduk.log

: <<'END'
# Remove duplicates - e.g. minidentity=95 allow 5 substitutions and 3 indels
dedupe.sh  -Xmx32g -eoom \
    in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
    in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
    out="${QC_FASTQ}/deduped-interleaved_trimmed.fastq.gz" \
    minidentity=95
    #s=5 e=5

reformat.sh -Xmx32g -eoom \
    in="${QC_FASTQ}/deduped-interleaved_trimmed.fastq.gz" \
    out1="${BAMQC_OUTPUT}/read1_filtered.fastq.gz" \
    out2="${BAMQC_OUTPUT}/read2_filtered.fastq.gz"

    ## Build the database
    # --threads ${coreNum} \
    bowtie2-build --seed 42 \
        ${ASSEMBLY_OUTPUT}/contigs.fasta \
        ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} |\
        tee -a ${LOG_DIR}/bowtie2_build_db_index.log

    ## Map the reads
    bowtie2 -t -D 10 -R 2 -L 31 -i S,0,2.50 -N 0 \
        -X ${maxInsert} \
        -k ${maxAlignments} \
        --threads ${coreNum} \
        -x ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
        --no-mixed \
        --no-discordant \
        --end-to-end \
        -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
        -2 "${QC_FASTQ}/read2_trimmed.fastq.gz" | \
    samtools view \
        -bh \
        -@ ${coreNum} \
        -f 0x003 \
        -o ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam - |\
        tee -a ${LOG_DIR}/bowtie2_mapping.log

    # Remove PCR duplicates
    samtools sort -n -o ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam
    samtools fixmate -cm ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.fixmate.bam
    samtools sort -o ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.fixmate.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.fixmate.bam
    samtools markdup -s -S ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.fixmate.bam ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.markdup.bam


END

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT_2} \
"${QC_FASTQ}/read1_trimmed.fastq.gz" \
"${QC_FASTQ}/read2_trimmed.fastq.gz"


## Run SPAdes assembly
timem spades.py \
-t ${coreNum} \
--careful \
-k 27,47,63,77,89,99,107,115,121,127 \
--pe1-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
--pe1-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} |\
tee -a ${LOG_DIR}/spades_assembly.log

#filter assembly contigs by Length - at least 1000
reformat.sh in="${ASSEMBLY_OUTPUT}/contigs.fasta" out="${ASSEMBLY_OUTPUT}/filtered_contigs_1kbp.fasta" minlength=1000

LOCAL_DB_NAME="contigs"
OUTPUT_PREFIX="trimmed_reads_vs_${LOCAL_DB_NAME}"

maxInsert=${maxInsert:-3000};
maxAlignments=${maxAlignments:-300};

## Build the database
# --threads ${coreNum} \
bowtie2-build --seed 42 \
    ${ASSEMBLY_OUTPUT}/contigs.fasta \
    ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} |\
    tee -a ${LOG_DIR}/bowtie2_build_db_index.log

bowtie2 --sensitive-local -p ${coreNum} --seed 42 -x ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" -2 "${QC_FASTQ}/read2_trimmed.fastq.gz" -S "${BWT_OUTPUT}/contigs_aligned.sam"
# index contigs
samtools faidx "${ASSEMBLY_OUTPUT}/contigs.fasta"
# create bam file
samtools import "${ASSEMBLY_OUTPUT}/contigs.fasta.fai" "${BWT_OUTPUT}/contigs_aligned.sam" "${BWT_OUTPUT}/contigs_aligned.bam"
# sort bam file
samtools sort -@ ${coreNum} "${BWT_OUTPUT}/contigs_aligned.bam"  -o "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# index bam
samtools index -@ ${coreNum}  "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# generate contigs stats
samtools idxstats "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.idxstats.txt"
samtools flagstat "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.flagstat.txt"

#bamqc
#qualimap bamqc -outdir "${BAMQC_OUTPUT}" -bam "${BWT_OUTPUT}/contigs_aligned.sorted.bam" -c
#qualimap bamqc -outdir "${BAMQC_OUTPUT}" -bam "${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.markdup.bam" -c

# Run Quast without reference
#--rna-finding \
quast.py \
-t ${coreNum} \
--glimmer \
--rna-finding \
--contig-thresholds 0,1000,5000,10000,25000,50000,100000,250000,500000,1000000 \
${ASSEMBLY_OUTPUT}/contigs.fasta \
-o ${QUAST_OUTPUT} | tee -a ${LOG_DIR}/quast.log

#--pe1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
#--pe2 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
