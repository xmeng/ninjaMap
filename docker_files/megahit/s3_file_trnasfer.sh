# Files transfer from s3 to sherlock
while IFS=$'\t' read -r -a line
do
SAMPLE=${line[0]}
echo $SAMPLE

rclone sync ${line[2]}/${line[1]}/ "${SCRATCH}/Samples/${SAMPLE}/"
#rm "${SCRATCH}/Samples/${SAMPLE}/*.zip" "${SCRATCH}/Samples/${SAMPLE}/*.html"
mv ${SCRATCH}/Samples/${SAMPLE}/*_R1.fastq.gz  ${SCRATCH}/Samples/${SAMPLE}/${SAMPLE}_R1_trimmed.fastq.gz
mv ${SCRATCH}/Samples/${SAMPLE}/*_R2.fastq.gz  ${SCRATCH}/Samples/${SAMPLE}/${SAMPLE}_R2_trimmed.fastq.gz

done < ~/handuo/test22.txt
