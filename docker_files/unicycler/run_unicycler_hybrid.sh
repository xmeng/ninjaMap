#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-16}
ShortTargetBase=${ShortTargetBase:-500000000}
LongTargetBase=${LongTargetBase:-500000000}

# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#longreads="${3}" input long reads in the sorted bam format
#S3OUTPUTPATH = "${4}"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/UNICYCLER"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc"
FASTQC_OUTPUT2="${LOCAL_OUTPUT}/fastqc_pacbio"
FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${FASTQ_NAME})

LOCAL_DB_PATH=${LOCAL}/databases
LOCAL_DB_NAME="contigs"
OUTPUT_PREFIX="trimmed_reads_vs_${LOCAL_DB_NAME}"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
BWT_OUTPUT="${LOCAL_OUTPUT}/bowtie2"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${FASTQC_OUTPUT2}"
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}"  "${BAMQC_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"
aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.bam"

#cp  ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#cp  ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

# Constant definitions for bbduk
adapterFile="adapters,phix"
trimQuality=${trimQuality:-25}
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim short reads, -eoom exits when out of memory
timem bbduk.sh -Xmx16g tbo -eoom hdist=1 qtrim=rl ktrim=r \
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

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT} \
"${QC_FASTQ}/read1_trimmed.fastq.gz" \
"${QC_FASTQ}/read2_trimmed.fastq.gz"

# convert bam to fastq for long reads
samtools bam2fq  "${RAW_FASTQ}/long.bam" >  "${RAW_FASTQ}/long.fastq"

# filter long reads with reference --target_bases 500000000 \
filtlong \
-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
--min_length 1000 \
--keep_percent 90 \
--target_bases ${LongTargetBase} \
--trim \
--split 500 \
"${RAW_FASTQ}/long.fastq" | gzip > "${QC_FASTQ}/long_trimmed.fastq.gz"


#Run fastqc for long reads
#fastqc \
#-t ${coreNum} \
#-o ${FASTQC_OUTPUT2} \
#"${QC_FASTQ}/long_trimmed.fastq.gz"

# Downsample the short reads if Uncycler is hang
reformat.sh \
samplebasestarget=${ShortTargetBase} \
in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
out="${QC_FASTQ}/read1_sampled.fastq.gz" \
out2="${QC_FASTQ}/read2_sampled.fastq.gz"



## Run assembly
timem unicycler \
-t ${coreNum} \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" \
-2 "${QC_FASTQ}/read2_sampled.fastq.gz" \
-l "${QC_FASTQ}/long_trimmed.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} |\
tee -a ${LOG_DIR}/unicycler_assembly.log.txt


## Build the database
# --threads ${coreNum} \
bowtie2-build --seed 42 \
    ${ASSEMBLY_OUTPUT}/assembly.fasta \
    ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} |\
    tee -a ${LOG_DIR}/bowtie2_build_db_index.log

bowtie2 --sensitive-local -p ${coreNum} --seed 42 -x ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" -2 "${QC_FASTQ}/read2_sampled.fastq.gz" -S "${BWT_OUTPUT}/contigs_aligned.sam"
# index contigs
samtools faidx "${ASSEMBLY_OUTPUT}/assembly.fasta"
# create bam file
samtools import "${ASSEMBLY_OUTPUT}/assembly.fasta.fai" "${BWT_OUTPUT}/contigs_aligned.sam" "${BWT_OUTPUT}/contigs_aligned.bam"
# sort bam file
samtools sort -@ ${coreNum} "${BWT_OUTPUT}/contigs_aligned.bam"  -o "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# index bam
samtools index -@ ${coreNum}  "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# generate contigs stats
samtools idxstats "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.idxstats.txt"
samtools flagstat "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.flagstat.txt"

#bamqc
qualimap bamqc -outdir "${BAMQC_OUTPUT}" -bam "${BWT_OUTPUT}/contigs_aligned.sorted.bam" -c


# Run Quast without reference
#-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
#-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
#--pacbio "${RAW_FASTQ}/long.fastq.gz" \
timem quast.py \
-t ${coreNum} \
--glimmer \
--rna-finding \
--contig-thresholds 1000,5000,10000,25000,50000,100000,250000,500000,1000000 \
${ASSEMBLY_OUTPUT}/assembly.fasta \
-o ${QUAST_OUTPUT} | tee -a ${LOG_DIR}/quast.log.txt


######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
