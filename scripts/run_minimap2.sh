#!/bin/bash -x

set -e
set -u
set -o pipefail

FASTQ="${1}"
DB="${2}"
THREADS=6
OUTPUTDIR="alignment"

FASTQ_NO_EXT=${FASTQ%.*}
DB_NO_EXT=${DB%%.*}
OUTPUT_PREFIX="${FASTQ_NO_EXT##*/}_vs_${DB_NO_EXT##*/}"

SAM_DIR="${OUTPUTDIR}/sam"
BAM_DIR="${OUTPUTDIR}/bam"

mkdir -p ${SAM_DIR} ${BAM_DIR}

minimap2 \
    -t ${THREADS} \
    -2 \
    -Q \
    --MD \
    --eqx \
    --sam-hit-only \
    -ax asm10 \
    -N 100 \
    -o "${SAM_DIR}/${OUTPUT_PREFIX}.sam" \
    "${DB}" \
    "${FASTQ}" &> "${OUTPUT_PREFIX}.aln.log"

samtools \
    view \
    -@ ${THREADS} \
    -bh "${SAM_DIR}/${OUTPUT_PREFIX}.sam" | \
    samtools sort \
    -@ ${THREADS} \
    -o "${BAM_DIR}/${OUTPUT_PREFIX}.sortedByCoord.bam"

rm "${SAM_DIR}/${OUTPUT_PREFIX}.sam"
