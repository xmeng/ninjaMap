# Setup base image
FROM continuumio/miniconda3:latest
# 4.8.2

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN conda install --freeze-installed -y nomkl numpy pandas \
      && conda install --freeze-installed -c conda-forge -y awscli biopython \
      && conda install --freeze-installed -c bioconda -y bowtie2 bedtools \
                              vcftools samtools==1.9 sambamba pysam \
                              pysamstats pybedtools bbmap \
      && conda clean -afy

RUN mkdir -p /mnt
WORKDIR /mnt

# Get Repo
COPY . .

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.8.2" \
      software.name="ninjaMap" \
      software.description="Strain abundance pipeline" \
      software.website="" \
      container.category="aligner"

RUN chmod -R +rx ./
