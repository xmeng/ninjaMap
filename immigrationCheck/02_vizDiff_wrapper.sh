#!/bin/bash -x
# shellcheck disable=SC2086

set -euo pipefail

ORG_NAME=${1%/}
COMPARE_WITH="${2}"
METHOD=${3:-"majority"}
METHOD=${METHOD,,}
S3_OUTPUT_BASE_PATH=${4:-"s3://czbiohub-microbiome/Sunit_Jain/scratch/immigrationCheck/output/v3_20200813_${METHOD}"}
S3_OUTPUT_BASE_PATH=${S3_OUTPUT_BASE_PATH%/}

# Create comparison list
bash 02a_create_comparison_list.sh ${ORG_NAME} ${COMPARE_WITH}

# Compare
python 02b_compute_strain_difference.py ${ORG_NAME} ${ORG_NAME}/compare.s3paths.list &> ${ORG_NAME}.02b.log

# Adjust 02b output location
OUTPUTDIR="${ORG_NAME}/${COMPARE_WITH}"
mkdir -p ${OUTPUTDIR}/prob_diff_arrays
mv ${ORG_NAME}*.npy ${OUTPUTDIR}/prob_diff_arrays/
mv ${ORG_NAME}.02b* ${OUTPUTDIR}/

# Input vs Invader calls
rm -rf ${ORG_NAME}/control_profiles.for_2d.list ${ORG_NAME}/samples_profiles.for_2d.list
for i in $(seq 1 5); do
    find ${OUTPUTDIR}/prob_diff_arrays/ -name "*${COMPARE_WITH}M${i}.prob_diff.npy" >> ${ORG_NAME}/control_profiles.for_2d.list
done

for i in $(seq 6 17); do
    find ${OUTPUTDIR}/prob_diff_arrays/ -name "*${COMPARE_WITH}M${i}.prob_diff.npy" >> ${ORG_NAME}/samples_profiles.for_2d.list
done

{
    find ${OUTPUTDIR}/prob_diff_arrays/ -name "*Pat1*.prob_diff.npy" | sort 
    find ${OUTPUTDIR}/prob_diff_arrays/ -name "*Pat2*.prob_diff.npy" | sort
    find ${OUTPUTDIR}/prob_diff_arrays/ -name "*Pat3*.prob_diff.npy" | sort 
} >> ${ORG_NAME}/samples_profiles.for_2d.list

find ${OUTPUTDIR}/prob_diff_arrays/ -name  "*Com*.prob_diff.npy" | sort >> ${ORG_NAME}/samples_profiles.for_2d.list

python 02d_compare_distributions_wViz.py \
    ${ORG_NAME} \
    ${METHOD} \
    ${ORG_NAME}/control_profiles.for_2d.list \
    ${ORG_NAME}/samples_profiles.for_2d.list &> ${ORG_NAME}.02d.log

# Move 02c outputs
mv ${ORG_NAME}.02* ${OUTPUTDIR}/

# Save output to S3
aws s3 sync --quiet ${ORG_NAME} \
    ${S3_OUTPUT_BASE_PATH}/${ORG_NAME}

# remove output dir
rm -rf ${ORG_NAME:?}/