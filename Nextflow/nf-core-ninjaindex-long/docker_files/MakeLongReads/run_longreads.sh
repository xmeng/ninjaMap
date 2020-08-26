#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

S3FASTA="${1}"
COV_FOLD=${2:-10}


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
PAIRED_FASTQ="${LOCAL_OUTPUT}/fastq"
LOCAL_FASTA="${OUTPUTDIR}/${FASTA}"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${PAIRED_FASTQ}"
#aws s3 cp ${S3FASTA} ${LOCAL_FASTA}
cp ${S3FASTA} ${LOCAL_FASTA}

LongReads="${PAIRED_FASTQ}/${PREFIX}.fastq.gz"



# Generate index first
bbmap.sh ref="${LOCAL_FASTA}"

# Generate Synthetic long reads command line
randomreads.sh \
in="${LOCAL_FASTA}" \
coverage=${COV_FOLD} \
maxsnps=0 \
insrate=0 \
adderrors=false \
paired=f \
out=${PREFIX}.perfect_long.fastq \
minlength=1000 \
maxlength=60000 \
seed=5

#compress fastq files
gzip < ${PREFIX}.perfect_long.fastq > $LongReads


# Sync
#aws s3 sync ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
#aws s3 cp ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
