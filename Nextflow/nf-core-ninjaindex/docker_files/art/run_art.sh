#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

S3FASTA="${1}"
#S3OUTPUTPATH="${2}"
#S3OUTPUTPATH=${S3OUTPUTPATH%/}
PIGZ_COMPRESSION_THREADS=${CORE_NUM:-4}

FASTA=$(basename -- "$S3FASTA")
PREFIX="${FASTA%.*}"

LOCAL=$(pwd)
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
PAIRED_FASTQ="${LOCAL_OUTPUT}/paired_fastq"
LOCAL_FASTA="${OUTPUTDIR}/${FASTA}"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${PAIRED_FASTQ}"
#aws s3 cp ${S3FASTA} ${LOCAL_FASTA}
cp ${S3FASTA} ${LOCAL_FASTA}

FWD="${PAIRED_FASTQ}/${PREFIX}.R1.fastq.gz"
REV="${PAIRED_FASTQ}/${PREFIX}.R2.fastq.gz"

# NUM_SEQS=$(grep -c ">" ${LOCAL_FASTA})
COV_FOLD=${COV_FOLD:-10}

# Grinder Command
/mnt/art_bin_MountRainier/art_illumina \
-na \
-ef \
-ss HS25 \
-i "${LOCAL_FASTA}" \
-p \
-l 150 \
-f ${COV_FOLD} \
-m 500 \
-s 10 \
-o "${PREFIX}" &> "${LOG_DIR}.log"

#compress fastq files
gzip < ${PREFIX}1.fq > $FWD
gzip < ${PREFIX}2.fq > $REV

# Sync
#aws s3 sync ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
#aws s3 cp ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
