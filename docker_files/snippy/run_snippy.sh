#!/bin/bash -x

set -euoE pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
DEFAULT=$coreNum
coreNum=${coreNum:-$DEFAULT}

mincov=${mincov:-10}
minfrac=${minfrac:-0}

if [ -z "${coreNum:-}" ]; then
  echo "coreNum was not set"
fi

# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#reference="${3}" reference genome
#insertion="${4}"   inserted gene
#S3OUTPUTPATH = "${5}"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSMEBLY_OUTPUT="${LOCAL_OUTPUT}/SPAdes"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT_1="${LOCAL_OUTPUT}/fastqc_raw"
FASTQC_OUTPUT_2="${LOCAL_OUTPUT}/fastqc_filtered"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${FASTQ_NAME})
BAMQC_OUTPUT2="${LOCAL_OUTPUT}/bamqc_reads_vs_reference"
MINIMAP2_OUTPUT="${LOCAL_OUTPUT}/bam_insertion_vs_contigs"
BAMQC_OUTPUT3="${LOCAL_OUTPUT}/bamqc_insertion_vs_contigs"
SNP_OUTPUT="${LOCAL_OUTPUT}/snippy"
DIR_NAME=${S3OUTPUTPATH##*/}

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${SNP_OUTPUT}"


hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"
aws s3 cp --quiet ${reference} "${RAW_FASTQ}/ref.fa"
aws s3 cp --quiet ${insertion} "${RAW_FASTQ}/insertion.fa"

# combine ref and inserted gene into one fasta file
cat "${RAW_FASTQ}/ref.fa" "${RAW_FASTQ}/insertion.fa" > "${RAW_FASTQ}/ref_combination.fa"

#cp  ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#cp  ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

: <<'END'
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
timem bbduk.sh -Xmx16g tbo -eoom hdist=1 qtrim=rl ktrim=r \
    entropy=0.5 entropywindow=50 entropyk=5 \
    in1="${RAW_FASTQ}/read1.fastq.gz" \
    in2="${RAW_FASTQ}/read2.fastq.gz" \
    out1="${BAMQC_OUTPUT}/read1_filtered.fastq.gz" \
    out2="${BAMQC_OUTPUT}/read2_filtered.fastq.gz" \
    ref=${adapterFile} \
    k="${kmer_value}" \
    mink="${min_kmer_value}" \
    trimq="${trimQuality}" \
    minlen="${minLength}" \
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" |\
    tee -a ${LOG_DIR}/bbduk.log


# Remove duplicates - e.g. minidentity=95 allow 5 substitutions and 3 indels
dedupe.sh  -Xmx16g -eoom \
    in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
    in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
    out="${QC_FASTQ}/deduped-interleaved_trimmed.fastq.gz" \
    minidentity=90
    #s=10 e=5

reformat.sh -Xmx16g -eoom \
    in="${QC_FASTQ}/deduped-interleaved_trimmed.fastq.gz" \
    out1="${BAMQC_OUTPUT}/read1_filtered.fastq.gz" \
    out2="${BAMQC_OUTPUT}/read2_filtered.fastq.gz"


#Align to reference
bwa mem -t ${coreNum}

# Remove PCR duplicates
samtools sort -n -o ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam
samtools fixmate -cm ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.fixmate.bam
samtools sort -o ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.fixmate.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.fixmate.bam
samtools markdup -s -r ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.fixmate.bam ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.dedupe.bam


samplereadstarget=${samplereadstarget:-0}
target=${target:-0}

#To sample 10% of the reads: samplerate=0.1
if [ $samplereadstarget -gt 0 ]
then
  reformat.sh -Xmx16g -eoom \
      in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
      in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
      out="${QC_FASTQ}/read1_trimmed_sampled.fastq.gz" \
      out2="${QC_FASTQ}/read2_trimmed_sampled.fastq.gz" \
      samplereadstarget="${samplereadstarget}"

#output file of reads with an average depth of 100x.
#Reads with an apparent depth of under 5x will be presumed to be errors and discarded.
elif [ $target -gt 0 ]
then
bbnorm.sh -Xmx16g -eoom \
  in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
  in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
  out="${QC_FASTQ}/read1_trimmed_normalized.fastq.gz" \
  out2="${QC_FASTQ}/read2_trimmed_normalized.fastq.gz" \
  target=$target \
  min=5
fi

if [ -e "${QC_FASTQ}/read1_trimmed_sampled.fastq.gz" ]
then
    ln -s "${QC_FASTQ}/read1_trimmed_sampled.fastq.gz" "${QC_FASTQ}/read1_filtered.fastq.gz"
    ln -s "${QC_FASTQ}/read1_trimmed_sampled.fastq.gz" "${QC_FASTQ}/read2_filtered.fastq.gz"

elif [ -e "${QC_FASTQ}/read1_trimmed_normalized.fastq.gz" ]
then
    ln -s "${QC_FASTQ}/read1_trimmed_normalized.fastq.gz" "${QC_FASTQ}/read1_filtered.fastq.gz"
    ln -s "${QC_FASTQ}/read2_trimmed_normalized.fastq.gz" "${QC_FASTQ}/read2_filtered.fastq.gz"
else
    ln -s "${QC_FASTQ}/read1_trimmed.fastq.gz"  "${QC_FASTQ}/read1_filtered.fastq.gz"
    ln -s "${QC_FASTQ}/read2_trimmed.fastq.gz"  "${QC_FASTQ}/read2_filtered.fastq.gz"
fi


#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT_2} \
"${BAMQC_OUTPUT}/read1_filtered.fastq.gz" \
"${BAMQC_OUTPUT}/read2_filtered.fastq.gz"


## Run SPAdes assembly 21,33,55,77,99,127
timem spades.py \
-t ${coreNum} \
--careful \
--untrusted-contigs "${RAW_FASTQ}/ref_combination.fa" \
-k 27,47,63,77,89,99,107,115,121,127 \
--pe1-1 "${BAMQC_OUTPUT}/read1_filtered.fastq.gz" \
--pe1-2 "${BAMQC_OUTPUT}/read2_filtered.fastq.gz" \
-o ${ASSMEBLY_OUTPUT} |\
tee -a ${LOG_DIR}/spades_assembly.log.txt

# Filter contig length by at leat 1000 bp
reformat.sh in="${ASSMEBLY_OUTPUT}/contigs.fasta" out="${ASSMEBLY_OUTPUT}/filtered_contigs_1kb.fasta" minlength=1000

#BWA alignment reads to contigs
bwa index "${ASSMEBLY_OUTPUT}/contigs.fasta"
bwa mem -t ${coreNum} "${ASSMEBLY_OUTPUT}/contigs.fasta" "${BAMQC_OUTPUT}/read1_filtered.fastq.gz" "${BAMQC_OUTPUT}/read2_filtered.fastq.gz" > "${BAMQC_OUTPUT}/reads_aligned_contigs.sam"
samtools view -bS "${BAMQC_OUTPUT}/reads_aligned_contigs.sam" > "${BAMQC_OUTPUT}/reads_aligned_contigs.bam"
samtools sort -o "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam" "${BAMQC_OUTPUT}/reads_aligned_contigs.bam"
samtools index "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam"
samtools flagstat "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam" > "${LOG_DIR}/reads_aligned_contigs.flagstat.txt"
# bam qc
qualimap bamqc -nt ${coreNum} -outdir "${BAMQC_OUTPUT}" -bam "${BAMQC_OUTPUT}/reads_aligned_contigs.sorted.bam" -c


#BWA alignment - reads to reference
bwa index "${RAW_FASTQ}/ref_combination.fa"
bwa mem -t ${coreNum} "${RAW_FASTQ}/ref_combination.fa"  "${BAMQC_OUTPUT}/read1_filtered.fastq.gz" "${BAMQC_OUTPUT}/read2_filtered.fastq.gz" > "${BAMQC_OUTPUT2}/reads_aligned_ref.sam"
samtools view -bS "${BAMQC_OUTPUT2}/reads_aligned_ref.sam" > "${BAMQC_OUTPUT2}/reads_aligned_ref.bam"
samtools sort -o "${BAMQC_OUTPUT2}/reads_aligned_ref.sorted.bam" "${BAMQC_OUTPUT2}/reads_aligned_ref.bam"
samtools index "${BAMQC_OUTPUT2}/reads_aligned_ref.sorted.bam"
samtools flagstat "${BAMQC_OUTPUT2}/reads_aligned_ref.sorted.bam" > "${LOG_DIR}/reads_aligned_ref.flagstat.txt"
# bam qc
qualimap bamqc -nt ${coreNum} -outdir "${BAMQC_OUTPUT2}" -bam "${BAMQC_OUTPUT2}/reads_aligned_ref.sorted.bam" -c

# align inserted gene to the contigs
cp "${ASSMEBLY_OUTPUT}/filtered_contigs_1kb.fasta"  ${MINIMAP2_OUTPUT}/${DIR_NAME}_contigs_1kb.fasta
cp  "${RAW_FASTQ}/insertion.fa"  ${MINIMAP2_OUTPUT}/${DIR_NAME}_insertion.fasta
minimap2 -t ${coreNum} -a "${ASSMEBLY_OUTPUT}/filtered_contigs_1kb.fasta"  "${RAW_FASTQ}/insertion.fa" | samtools sort -@${coreNum} -o "${MINIMAP2_OUTPUT}/${DIR_NAME}_insertion_vs_contigs.sorted.bam" -
samtools index "${MINIMAP2_OUTPUT}/${DIR_NAME}_insertion_vs_contigs.sorted.bam"
# bam qc
qualimap bamqc -nt ${coreNum} -outdir "${BAMQC_OUTPUT3}" -bam "${MINIMAP2_OUTPUT}/${DIR_NAME}_insertion_vs_contigs.sorted.bam" -c

END

#Rapid haploid variant calling and core genome alignment
snippy --cpus ${coreNum} \
--report \
--force \
--mincov ${mincov} \
--minfrac ${minfrac} \
--outdir "${SNP_OUTPUT}" \
--ref "${RAW_FASTQ}/ref_combination.fa" \
--R1 "${RAW_FASTQ}/read1.fastq.gz" \
--R2 "${RAW_FASTQ}/read2.fastq.gz" | tee -a ${LOG_DIR}/snippy.log.txt


######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
