#!/bin/bash -x

set -euoE pipefail

ulimit -c unlimited

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
#longreads="${3}" input long reads
#S3OUTPUTPATH = "${4}"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"

DATABASE="${LOCAL}/db"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
QC_FASTQ="${LOCAL_OUTPUT}/trimmed_fastq"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSMEBLY_OUTPUT="${LOCAL_OUTPUT}/metaSPAdes"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc"
FASTQC_OUTPUT2="${LOCAL_OUTPUT}/fastqc2"
FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${FASTQ_NAME})

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${ASSMEBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${DATABASE}" "${FASTQC_OUTPUT2}"
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



# Remove host contamination - human
# 1. Download human genome ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38//GCA_000001405.15_GRCh38_genomic.fna.gz
# rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz .
#curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38//GCA_000001405.15_GRCh38_genomic.fna.gz
#gunzip GCA_000001405.15_GRCh38_genomic.fna.gz
#mv GCA_000001405.15_GRCh38_genomic.fna "${DATABASE}/human.fna"
# 2. build human genome database
#bowtie2-build --large-index "${DATABASE}/human.fna" "${DATABASE}/human" --seed 42 --verbose | tee -a ${LOG_DIR}/bowtie2_build_db_human_index.log
# copy human index from s3
aws s3 sync s3://czbiohub-microbiome/Xiandong_Meng/database/hg38/ "${DATABASE}"
# 3. map read to human database
bowtie2 --very-sensitive-local -p  ${coreNum} --seed 42 -x "${DATABASE}/human" -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" -2 "${QC_FASTQ}/read1_trimmed.fastq.gz" -S "${QC_FASTQ}/human.sam"
# 4. convert sam to bam
samtools view -@ ${coreNum} -bS "${QC_FASTQ}/human.sam" > "${QC_FASTQ}/human.bam"
# 5. Extract unmapped reads
#    -f 12 Extracts only alignments with both reads unmapped
#    -F 256 Does not extract aligns that are not primary alignment
samtools view -@ ${coreNum} -b -f 12 -F 256 "${QC_FASTQ}/human.bam" > "${QC_FASTQ}/unmapped_human.bam"
# 6. Sort BAM file
samtools sort -n "${QC_FASTQ}/unmapped_human.bam" -o "${QC_FASTQ}/unmapped_human.sorted.bam"
# 7. convert unmapped reads in bam to fastq
samtools bam2fq "${QC_FASTQ}/unmapped_human.sorted.bam" > "${QC_FASTQ}/cleaned_no_human.fastq"

: <<'END'

# Remove host contamination - mouse
# 1. Download mouse genome
#curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_GRCm38.p6/GCA_000001635.8_GRCm38.p6_genomic.fna.gz
#gunzip GCA_000001635.8_GRCm38.p6_genomic.fna.gz
#mv GCA_000001635.8_GRCm38.p6_genomic.fna "${DATABASE}/mouse.fna"
## 2. build mouse genome database
#bowtie2-build --large-index "${DATABASE}/mouse.fna" "${DATABASE}/mouse" --seed 42 --verbose | tee -a ${LOG_DIR}/bowtie2_build_db_mouse_index.log
# copy human index from s3
aws s3 sync s3://czbiohub-microbiome/Xiandong_Meng/database/mouse/ "${DATABASE}"
# 3. map interleaved read to mouse database
bowtie2 --very-sensitive-local -p ${coreNum} --seed 42 -x "${DATABASE}/mouse" --12 "${QC_FASTQ}/cleaned_no_human.fastq" -S "${QC_FASTQ}/mouse.sam"
# 4. convert sam to bam
samtools view -@ ${coreNum} -bS "${QC_FASTQ}/mouse.sam" > "${QC_FASTQ}/mouse.bam"
# 5. Extract unmapped reads
#    -f 12 Extracts only alignments with both reads unmapped
#    -F 256 Does not extract aligns that are not primary alignment
samtools view -@ ${coreNum} -b -f 12 -F 256 "${QC_FASTQ}/mouse.bam" > "${QC_FASTQ}/unmapped_mouse.bam"
# 6. Sort BAM file
samtools sort -n "${QC_FASTQ}/unmapped_mouse.bam" -o "${QC_FASTQ}/unmapped_mouse.sorted.bam"
# 7. convert unmapped reads in bam to fastq
samtools bam2fq "${QC_FASTQ}/unmapped_mouse.sorted.bam" > "${QC_FASTQ}/cleaned_no_human_mouse.fastq"

END

#Run fastqc for cleaned reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT2} \
"${QC_FASTQ}/cleaned_no_human.fastq"
#"${QC_FASTQ}/cleaned_no_human_mouse.fastq"




######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
