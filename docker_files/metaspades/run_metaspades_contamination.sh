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
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
DATABASE="${LOCAL}/db"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
QC_FASTQ="${LOCAL_OUTPUT}/trimmed_fastq"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/metaSPAdes"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc-filtering"
FASTQC_OUTPUT2="${LOCAL_OUTPUT}/fastqc-decontamination"
FASTQC_OUTPUT3="${LOCAL_OUTPUT}/fastqc-deduplication"
FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${FASTQ_NAME})

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${DATABASE}" "${FASTQC_OUTPUT2}" "${FASTQC_OUTPUT3}"
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

#reformat.sh \
#-in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
#-in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
#-out="${QC_FASTQ}/interleaved_trimmed.fastq.gz"

# Remove duplicates - e.g. minidentity=95 allow 5 substitutions and 3 indels
dedupe.sh  -Xmx16g -eoom \
in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
out="${QC_FASTQ}/deduped-interleaved_trimmed.fastq.gz" \
s=5 e=3

#Run fastqc for deduped short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT3} \
"${QC_FASTQ}/deduped-interleaved_trimmed.fastq.gz"



#: <<'END'
# Remove host contamination - human
# 1. Download human genome ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38//GCA_000001405.15_GRCh38_genomic.fna.gz
# rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz .
#curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38//GCA_000001405.15_GRCh38_genomic.fna.gz
#gunzip GCA_000001405.15_GRCh38_genomic.fna.gz
#mv GCA_000001405.15_GRCh38_genomic.fna "${DATABASE}/human.fna"
# 2. build human genome database
#bowtie2-build --large-index "${DATABASE}/human.fna" "${DATABASE}/human" --seed 42 --verbose | tee -a ${LOG_DIR}/bowtie2_build_db_human_index.log
# copy human index from s3
aws s3 sync s3://czbiohub-microbiome/Xiandong_Meng/database/BWA/ "${DATABASE}"
# 3. map read to human database --very-sensitive-local
#bowtie2  -p  ${coreNum} --seed 42 -x "${DATABASE}/human" -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" -2 "${QC_FASTQ}/read1_trimmed.fastq.gz" -S "${QC_FASTQ}/human.sam"
bwa mem -t ${coreNum}  "${DATABASE}/human.fna" -p "${QC_FASTQ}/deduped-interleaved_trimmed.fastq.gz" > "${QC_FASTQ}/bwa-human.sam"
#bwa mem -t ${coreNum}  "${DATABASE}/human.fna" "${QC_FASTQ}/read1_trimmed.fastq.gz" "${QC_FASTQ}/read1_trimmed.fastq.gz" > "${QC_FASTQ}/bwa-human.sam"
# 4. convert sam to bam
#samtools view -@ ${coreNum} -bS "${QC_FASTQ}/human.sam" > "${QC_FASTQ}/human.bam"
# 5. Extract unmapped reads
#    -f 12 Extracts only alignments with both reads unmapped
#    -F 256 Does not extract aligns that are not primary alignment
#samtools view -@ ${coreNum} -b -f 12 -F 256 "${QC_FASTQ}/human.bam" > "${QC_FASTQ}/unmapped_human.bam"
# 6. Sort, index and stat of BAM file
#samtools sort -n "${QC_FASTQ}/unmapped_human.bam" -o "${QC_FASTQ}/unmapped_human.sorted.bam"
#samtools index "${QC_FASTQ}/unmapped_human.sorted.bam"

samtools view -Sh -f4 "${QC_FASTQ}/bwa-human.sam" > "${QC_FASTQ}/bwa-human.unmapped.sam"
samtools flagstat "${QC_FASTQ}/bwa-human.unmapped.sam" > "${LOG_DIR}/bwa-human_mapping_stats.txt"

# 7. convert unmapped reads in bam to fastq
samtools bam2fq "${QC_FASTQ}/bwa-human.unmapped.sam" > "${QC_FASTQ}/cleaned_no_human.fastq"
# drop unpaired from interleaved PE fastq
seqtk dropse "${QC_FASTQ}/cleaned_no_human.fastq"  > "${QC_FASTQ}/cleaned_no_human.pe.fastq"


# Remove host contamination - mouse
# 1. Download mouse genome
#curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_GRCm38.p6/GCA_000001635.8_GRCm38.p6_genomic.fna.gz
#gunzip GCA_000001635.8_GRCm38.p6_genomic.fna.gz
#mv GCA_000001635.8_GRCm38.p6_genomic.fna "${DATABASE}/mouse.fna"
## 2. build mouse genome database
#bowtie2-build --large-index "${DATABASE}/mouse.fna" "${DATABASE}/mouse" --seed 42 --verbose | tee -a ${LOG_DIR}/bowtie2_build_db_mouse_index.log
# copy human index from s3
#aws s3 sync s3://czbiohub-microbiome/Xiandong_Meng/database/bwa/mouse/ "${DATABASE}"
# 3. map interleaved read to mouse database --very-sensitive-local
#bowtie2  -p ${coreNum} --seed 42 -x "${DATABASE}/mouse" --12 "${QC_FASTQ}/cleaned_no_human.fastq" -S "${QC_FASTQ}/mouse.sam"
bwa mem -t ${coreNum}  "${DATABASE}/mouse.fna" -p "${QC_FASTQ}/cleaned_no_human.pe.fastq" > "${QC_FASTQ}/bwa-mouse.sam"
# 4. convert sam to bam
#samtools view -@ ${coreNum} -bS "${QC_FASTQ}/mouse.sam" > "${QC_FASTQ}/mouse.bam"
# 5. Extract unmapped reads
#    -f 12 Extracts only alignments with both reads unmapped
#    -F 256 Does not extract aligns that are not primary alignment
#samtools view -@ ${coreNum} -b -f 12 -F 256 "${QC_FASTQ}/mouse.bam" > "${QC_FASTQ}/unmapped_mouse.bam"
# 6. Sort, index and stat of BAM file
#samtools sort -n "${QC_FASTQ}/unmapped_mouse.bam" -o "${QC_FASTQ}/unmapped_mouse.sorted.bam"
#samtools index "${QC_FASTQ}/unmapped_mouse.sorted.bam"

samtools view -Sh -f4 "${QC_FASTQ}/bwa-mouse.sam" > "${QC_FASTQ}/bwa-mouse.unmapped.sam"
samtools flagstat "${QC_FASTQ}/bwa-mouse.unmapped.sam" > "${LOG_DIR}/bwa-mouse_mapping_stats.txt"

# 7. convert unmapped reads in bam to fastq
samtools bam2fq "${QC_FASTQ}/bwa-mouse.unmapped.sam" > "${QC_FASTQ}/cleaned_no_human_mouse.fastq"
# drop unpaired from interleaved PE fastq
seqtk dropse "${QC_FASTQ}/cleaned_no_human_mouse.fastq"  > "${QC_FASTQ}/cleaned_no_human_mouse.pe.fastq"

#Run fastqc for cleaned reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT2} \
"${QC_FASTQ}/cleaned_no_human_mouse.pe.fastq"

#END

# convert bam to fastq for long reads
samtools bam2fq  "${RAW_FASTQ}/long.bam" >  "${RAW_FASTQ}/long.fastq"

# filter long reads #--target_bases 500000000 \
filtlong \
-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
--min_length 1000 \
--keep_percent 90 \
--trim \
--split 500 \
"${RAW_FASTQ}/long.fastq" | gzip > "${QC_FASTQ}/long_trimmed.fastq.gz"

#if using the pacbio clr reads
#--pacbio pacbio_clr.fastq

## Run metaSPAdes assembly Don't work with --trusted-contigs or --untrusted-contigs option
timem spades.py --meta \
-t ${coreNum} \
-k 21,33,55,77,99,127 \
--12 "${QC_FASTQ}/cleaned_no_human_mouse.pe.fastq" \
--pacbio "${QC_FASTQ}/long_trimmed.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} |\
tee -a ${LOG_DIR}/metaspades_assembly.log


# Run Quast without reference
#/bin/bash -c "source activate quast" &&
metaquast.py \
-t ${coreNum} \
--12 "${QC_FASTQ}/cleaned_no_human_mouse.pe.fastq" \
--pacbio "${QC_FASTQ}/long_trimmed.fastq.gz" \
${ASSEMBLY_OUTPUT}/contigs.fasta \
-o ${QUAST_OUTPUT} | tee -a ${LOG_DIR}/quast.log

#--glimmer \
#--rna-finding \
#-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
#-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
