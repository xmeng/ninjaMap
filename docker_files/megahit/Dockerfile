# Base Image
FROM continuumio/miniconda3:4.6.14
LABEL authors="Xiandong Meng" \
      description="Docker image for metagenome assembly " \
      software.name="megahit_mapping_binning" \
      software.description="An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph." \
      software.website="https://github.com/voutcn/megahit" \
      container.category="tool"

# Install system-level programs for quast
#RUN apt-get update && apt-get install -y zlib1g-dev
RUN apt-get update && apt-get install -y libopenblas-base pkg-config libfreetype6-dev libpng-dev python-matplotlib

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install software packages
RUN conda install --freeze-installed -c conda-forge -y awscli \
  && conda install --freeze-installed -c conda-forge -c bioconda -y megahit=1.2.9 python=3 bbmap samtools bowtie2 metabat2 libdeflate fastqc quast=5.0.2 \
  && conda clean -afy

# Install time and memory usage tool (timem)   qualimap
RUN conda install -c jrmadsen -y timemory


ADD run_megahit_mapping_binning.sh /usr/local/bin/run_megahit_mapping_binning.sh
RUN chmod +x /usr/local/bin/run_megahit_mapping_binning.sh

ADD run_bowtie2_mapping.sh /usr/local/bin/run_run_bowtie2_mapping.sh
RUN chmod +x /usr/local/bin/run_run_bowtie2_mapping.sh


ADD get_count_table.py /usr/local/bin/get_count_table.py
RUN chmod +x /usr/local/bin/get_count_table.py

RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
