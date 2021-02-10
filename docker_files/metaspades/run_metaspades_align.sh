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
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/metaSPAdes"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc"
FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${FASTQ_NAME})

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

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



#if using the pacbio clr reads
#--pacbio pacbio_clr.fastq

## Run metaSPAdes assembly Don't work with --trusted-contigs or --untrusted-contigs option
timem spades.py --meta \
-t ${coreNum} \
-k 21,33,55,77,99,127 \
--pe1-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
--pe1-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} |\
tee -a ${LOG_DIR}/metaspades_assembly.log

# Run Quast without reference
#/bin/bash -c "source activate quast" && metaquast.py \
metaquast \
-t ${coreNum} \
--glimmer \
--rna-finding \
--contig-thresholds 0,1000,5000,10000,25000,50000,100000,250000 \
-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
${ASSEMBLY_OUTPUT}/contigs.fasta \
-o ${QUAST_OUTPUT} | tee -a ${LOG_DIR}/mataquast.log


#BWA alignment reads to contigs
bwa index "${ASSEMBLY_OUTPUT}/contigs.fasta"
bwa mem  -t ${coreNum} "${ASSEMBLY_OUTPUT}/contigs.fasta" "${QC_FASTQ}/read1_trimmed.fastq.gz" "${QC_FASTQ}/read2_trimmed.fastq.gz" > "${BAMQC_OUTPUT}/reads_aligned_contigs.sam"
samtools view -bS "${BAMQC_OUTPUT}/reads_aligned_contigs.sam" > "${BAMQC_OUTPUT}/reads_aligned_contigs.bam"
samtools sort -o "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam" "${BAMQC_OUTPUT}/reads_aligned_contigs.bam"
samtools index "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam"
samtools flagstat "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam" > "${LOG_DIR}/reads_aligned_contigs.flagstat.txt"
# bam qc
qualimap bamqc -nt ${coreNum} -outdir "${BAMQC_OUTPUT}" -bam "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam" -c

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
