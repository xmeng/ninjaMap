#########################################
##     LAST UPDATED ON: 2020-04-24     ##
## NOT INTENDED TO BE USED AS A SCRIPT ##
##     SIMPLY A LOG OF THE WORKFLOW    ##
#########################################

# Create a list of genomes with their S3/FTP paths.
aws s3 ls s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422/source/original_genomes/ | awk -v s3path="s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422/fasta/original_genomes" '{printf "%s\t%s/%s\n", $4,s3path,$4}'| sed "s/.fasta.gz//" > original_genomes.s3paths.list 
aws s3 ls s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422/source/updated_genomes/ | awk -v s3path="s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422/fasta/updated_genomes" '{printf "%s\t%s/%s\n", $4,s3path,$4}'| sed "s/.fasta.gz//" > updated_genomes.s3paths.list 
cat original_genomes.s3paths.list updated_genomes.s3paths.list | wc -l
# 101
echo -e "Strain_Name\tFTP_link" > scv1_1.s3paths.list
cat original_genomes.s3paths.list updated_genomes.s3paths.list >> scv1_1.s3paths.list

# Standardize the genome names, sequence headers, generate metadata and concatenate each genome into a single file.
# python create_ninjamap_db.py genome_paths.list db_dir db_name
python create_ninjamap_db.py -list scv1_1.s3paths.list -db SCv1_1 &>create_ninjamap_db.log

aws s3 sync SCv1_1 s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422

cd SCv1_1 || exit 1

# Run Genome QC
bash -x run_binqc.sh qc fasta fna &> run_binqc.log

# Keep Genes and Proteins
mkdir -p qc/orfs/{genes,proteins}
find qc/GTDBtk/ -name '*_protein.fna' | parallel "cp {} qc/orfs/genes/"
find qc/GTDBtk/ -name '*_protein.fna' | parallel "cp {} qc/orfs/proteins/"

## Optionally, delete all GTDBtk intermediate results
# sudo chmod a+wx qc/GTDBtk/*/intermediate_results
# sudo rm -rf qc/GTDBtk/*/intermediate_results

grep -c '>' qc/orfs/genes/*_protein.fna | \
    sed -e 's/_protein.fna:/\t/' -e "s#qc/orfs/genes/##" > qc/orfs/num_genes.txt

aws s3 sync . s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422

cd db || exit 1

# Create a bowtie2 db of the concatenated fasta
docker container run --rm \
    --workdir $(pwd) \
    --volume $(pwd):$(pwd) \
    quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0 \
        bowtie2-build \
        --threads 8 \
        --seed 1712 \
        -f SCv1_1.fna \
        SCv1_1

# Cleanup
mkdir -p bowtie2_index
mv *bt2 bowtie2_index/

# Sync to S3
aws s3 sync . s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422/db

cd ../

python generate_ninjamap_db_report.py . SCv1.1
aws s3 sync . s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422

cd ../
