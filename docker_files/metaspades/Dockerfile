# Base Image
FROM continuumio/miniconda3:4.6.14
LABEL authors="Xiandong Meng" \
      description="Docker image for Hybrid assembly " \
      software.name="metaspades" \
      software.description="SPAdes – St. Petersburg genome assembler – is an assembly toolkit containing various assembly pipelines." \
      software.website="http://cab.spbu.ru" \
      container.category="tool"

# Install system-level programs for quast
#RUN apt-get update && apt-get install -y zlib1g-dev
RUN apt-get update && apt-get install -y libopenblas-base

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install software packages snippy
RUN conda install -c bioconda -y spades bbmap samtools fastqc bowtie2 bwa quast qualimap seqtk python=3 minimap2
#RUN conda install -c conda-forge -c bioconda -c defaults -y -c bioconda snippy bcftools
RUN conda install -c conda-forge -y awscli

# install long reads filter
RUN conda install -c bioconda filtlong

# Install time and memory usage tool (timem)
RUN conda install -c jrmadsen -y timemory

#COPY environment.yaml ./
#RUN conda env create -f environment.yaml && conda clean -a

#RUN /bin/bash -c "source activate quast"
#ENV PATH /opt/conda/envs/quast/bin:$PATH

# https://github.com/ablab/quast/issues/8
COPY search_references_meta.py /opt/conda/lib/python3.6/site-packages/quast_libs/search_references_meta.py

ADD run_metaspades.sh /usr/local/bin/run_metaspades.sh
RUN chmod +x /usr/local/bin/run_metaspades.sh

ADD run_metaspades_vital.sh /usr/local/bin/run_metaspades_vital.sh
RUN chmod +x /usr/local/bin/run_metaspades_vital.sh

ADD run_metaspades_short.sh /usr/local/bin/run_metaspades_short.sh
RUN chmod +x /usr/local/bin/run_metaspades_short.sh

ADD run_remove_contamination.sh /usr/local/bin/run_remove_contamination.sh
RUN chmod +x /usr/local/bin/run_remove_contamination.sh

ADD run_metaspades_contamination.sh /usr/local/bin/run_metaspades_contamination.sh
RUN chmod +x /usr/local/bin/run_metaspades_contamination.sh

ADD run_metaspades_align.sh /usr/local/bin/run_metaspades_align.sh
RUN chmod +x /usr/local/bin/run_metaspades_align.sh

ADD run_spades_ref.sh /usr/local/bin/run_spades_ref.sh
RUN chmod +x /usr/local/bin/run_spades_ref.sh

ADD run_spades_denovo.sh /usr/local/bin/run_spades_denovo.sh
RUN chmod +x /usr/local/bin/run_spades_denovo.sh

ADD run_spades_gene_insertion.sh /usr/local/bin/run_spades_gene_insertion.sh
RUN chmod +x /usr/local/bin/run_spades_gene_insertion.sh

ADD run_spades_hybrid_nanopore.sh /usr/local/bin/run_spades_hybrid_nanopore.sh
RUN chmod +x /usr/local/bin/run_spades_hybrid_nanopore.sh

RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
