#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-15};
LOCAL_DB_PATH=${LOCAL}/databases

S3DBPATH="${1}"
fastq="${2}"
#S3OUTPUTPATH

#S3DBPATH=s3://czbiohub-microbiome/Synthetic_Community/Genome_References/ncbi_fasta/Dorea-longicatena-DSM-13814-GCF_000154065.1_ASM15406v1.fna
#fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
#fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
#S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/minimap2_Test/Dorea-longicatena-DSM-13814

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
TMP_BWT_OUTPUT="${OUTPUTDIR}/minimap2"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
BWT_OUTPUT="${LOCAL_OUTPUT}/minimap2"
#S3OUTPUTPATH=${S3OUTPUTPATH%/}
S3DBPATH=${S3DBPATH%/*}
SAMPLE_NAME=$(basename ${fastq})

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" 
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}" "${TMP_BWT_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
#aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

LOCAL_REFSEQ_NAME=$(basename ${S3DBPATH})
LOCAL_REFSEQ_EXT=${S3DBPATH##*.}
LOCAL_DB_NAME=$(basename ${S3DBPATH} .${LOCAL_REFSEQ_EXT})

cp ${S3DBPATH} ${LOCAL_DB_PATH}/
cp ${fastq} "${RAW_FASTQ}/longread.fastq.gz"


OUTPUT_PREFIX="${SAMPLE_NAME}_vs_db_${LOCAL_DB_NAME}"

## long sequences against a reference genome
minimap2 \
    -t ${coreNum} \
    -a \
    ${LOCAL_DB_PATH}/${LOCAL_REFSEQ_NAME} \
    "${RAW_FASTQ}/longread.fastq.gz" | \
samtools view \
    -bh \
    -@ ${coreNum} \
    -o ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam - | \
    tee -a ${LOG_DIR}/read_mapping.log


# Remove PCR duplicates
samtools sort -n -o ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam
samtools markdup -s -S ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam  ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.markdup.bam
samtools index ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.markdup.bam

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
#aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
