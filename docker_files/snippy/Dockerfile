# Base Image
#FROM nfcore/base
FROM continuumio/miniconda3:4.6.14

LABEL authors="Xiandong Meng" \
      description="Docker image for snippy" \
      software.name="snippy" \
      software.description="Rapid haploid variant calling and core genome alignment" \
      software.website="https://github.com/tseemann/snippy" \
      container.category="tool"

# Update conda to latest version. /opt/conda/bin/
RUN conda update -n base -c defaults conda

# Install software packages
RUN conda install -c conda-forge -c defaults -y -c bioconda snippy

RUN conda install -c conda-forge -y awscli

# Install time and memory usage tool (timem)
#RUN /opt/conda/bin/conda install -c jrmadsen -y timemory

ADD run_snippy.sh /usr/bin/run_snippy.sh
RUN chmod +x /usr/bin/run_snippy.sh


RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
