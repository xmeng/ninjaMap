#!/bin/bash -x

set -euoE pipefail

START_TIME=$SECONDS
#export PATH="/opt/conda/bin:${PATH}"

LOCAL=/scratch/users/xdmeng/workdir
DEFAULT=20
coreNum=${coreNum:-$DEFAULT}

if [ -z "${coreNum:-}" ]; then
  echo "coreNum was not set"
fi

# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"

#group="${1}"
#sample="${2}"
#MYOUTPUTPATH="${3}"

MYOUTPUTPATH=/scratch/users/xdmeng/Analysis/${experiment}/comapping/${group}__${sample}
DB=/scratch/users/xdmeng/Analysis/${experiment}/${group}/bowtie2_index
FASTQ1=/scratch/users/xdmeng/Samples/${experiment}/${sample}/${sample}_R1_trimmed.fastq.gz
FASTQ2=/scratch/users/xdmeng/Samples/${experiment}/${sample}/${sample}_R2_trimmed.fastq.gz
# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_${1} #$( date +"%Y%m%d_%H%M%S" )
#RAW_FASTQ="${OUTPUTDIR}/raw_fastq"

LOCAL_TMP="${OUTPUTDIR}/tmp"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
MAPPING_OUTPUT="${LOCAL_OUTPUT}/mapping"
#FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc"
#BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_mapping"
#BINNING_OUTPUT="${LOCAL_OUTPUT}/bins"
#FASTQ_NAME=${fastq1%/*}
#SAMPLE_NAME=$(basename ${FASTQ_NAME})


mkdir -p ${OUTPUTDIR} ${LOCAL_OUTPUT} ${LOG_DIR}
mkdir -p ${MAPPING_OUTPUT} ${LOCAL_TMP}  ${MYOUTPUTPATH}
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

if [ -f $MYOUTPUTPATH/mapping/assembly_aligned.sorted.bam ]; then
    exit 0
fi

: <<'END'
for ((idx=0; idx<${#sampleList[@]}; ++idx));
do
    SAMPLE="${sampleList[idx]}"

    fastq1="s3://czbiohub-microbiome/Huang_Lab/Katharine_Ng/190906_Mouse_Metagenome/${SAMPLE}/03_NODUP/NODUP_HMN_UNMAPPED_TRIM_MARKED_${SAMPLE}_R1.fastq.gz"
    fastq2="s3://czbiohub-microbiome/Huang_Lab/Katharine_Ng/190906_Mouse_Metagenome/${SAMPLE}/03_NODUP/NODUP_HMN_UNMAPPED_TRIM_MARKED_${SAMPLE}_R2.fastq.gz"
    # Copy fastq.gz files from S3, only 2 files per sample
    aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/${SAMPLE}_R1.fastq.gz"
    aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/${SAMPLE}_R2.fastq.gz"
done

cat "${RAW_FASTQ}/*_R1.fastq.gz" > "${RAW_FASTQ}/read1.fastq.gz"
cat "${RAW_FASTQ}/*_R2.fastq.gz" > "${RAW_FASTQ}/read2.fastq.gz"


hash_kmer=${hash_kmer:-51}

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
    tossbrokenreads=t \
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" |\
    tee -a ${LOG_DIR}/bbduk.log

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT} \
"${QC_FASTQ}/read1_trimmed.fastq.gz" \
"${QC_FASTQ}/read2_trimmed.fastq.gz"

# downsample 80% of trimmed reads
reformat.sh \
samplerate=${samplerate} \
in="${QC_FASTQ}/read1_trimmed.fastq.gz" in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
out="${QC_FASTQ}/read1_sampled.fastq.gz" out2="${QC_FASTQ}/read2_sampled.fastq.gz" |\
tee -a ${LOG_DIR}/reformat.log.txt

END

# 2. Run mapping for short reads
# build bowtie index
#bowtie2-build  "${ASSEMBLY_OUTPUT}/filtered_contigs_1kb.fasta" "${DB}/database" –-seed 42

# map reads to database
bowtie2 --sensitive-local -p ${coreNum} --seed 42 -x "${DB}/database" \
-1 ${FASTQ1} -2 ${FASTQ2} | \
samtools view -@ ${coreNum} -bh -o "${LOCAL_TMP}/assembly_aligned.bam" - | \
tee -a ${LOG_DIR}/bowtie2_mapping.log.txt

# sort bam file
samtools sort -@ ${coreNum} "${LOCAL_TMP}/assembly_aligned.bam"  -o "${MAPPING_OUTPUT}/assembly_aligned.sorted.bam"
# index bam
samtools index -@ ${coreNum} "${MAPPING_OUTPUT}/assembly_aligned.sorted.bam"
# generate assembly stats
samtools idxstats "${MAPPING_OUTPUT}/assembly_aligned.sorted.bam" > "${LOG_DIR}/assembly_aligned.idxstats.txt"
samtools flagstat "${MAPPING_OUTPUT}/assembly_aligned.sorted.bam" > "${LOG_DIR}/assembly_aligned.flagstat.txt"
# check reads coverage/depth on contigs
#pileup.sh in="${MAPPING_OUTPUT}/assembly_aligned.sorted.bam" covstats="${LOG_DIR}/assembly_aligned.cov_stats.txt" > "${LOG_DIR}/assembly_aligned.stats.txt"

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
mv ${LOCAL_OUTPUT}/* ${MYOUTPUTPATH}
rm -rf ${OUTPUTDIR}
