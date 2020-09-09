#!/bin/bash -x
# shellcheck disable=SC2034
# shellcheck disable=SC2086
# shellcheck disable=SC2016

set -e
set -u
set -o pipefail

PREFIX=$1
REF=$2
QUERY=$3

DOCKER_IMAGE="quay.io/biocontainers/mummer"
DOCKER_VERSION="3.23--pl526he1b5a44_11"

DOCKER_CMD_PREFIX='docker container run --rm --workdir "$(pwd)" --volume "$(pwd)":"$(pwd)" ${DOCKER_IMAGE}:${DOCKER_VERSION}'

eval ${DOCKER_CMD_PREFIX} \
        nucmer \
        --maxgap=500 \
        --mincluster=100 \
        --prefix=${PREFIX} \
        ${REF} \
        ${QUERY} &> ${PREFIX}.nucmer.log

eval ${DOCKER_CMD_PREFIX} \
    show-coords \
    -r ${PREFIX}.delta > ${PREFIX}.coords

eval ${DOCKER_CMD_PREFIX} \
    delta-filter \
    -q \
    -r ${PREFIX}.delta > ${PREFIX}.filter

eval ${DOCKER_CMD_PREFIX} \
    mummerplot \
    ${PREFIX}.filter \
    -R ${REF} \
    -Q ${QUERY} \
    -t png \
    -s large \
    -p ${PREFIX}

# yum install gnuplot
gnuplot \
    -c ${PREFIX}.gp