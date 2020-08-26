FROM nfcore/base

LABEL authors="Xiandong Meng" \
      description="Docker image containing all requirements for nf-core/ninjaindex-long pipeline" \
      software="art" \
      about.summary="Generate synthetic next-generation sequencing reads." \
      about.home="bbtools" \
      software.version="bbmap-38.73-0"


USER root:root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y awscli pigz vim htop && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/*

RUN conda install -c bioconda -y bbmap
RUN conda clean -ya

RUN mkdir -p /mnt
WORKDIR /mnt


#RUN curl https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz --output - | tar xz
WORKDIR /mnt

ADD run_longreads.sh /usr/local/bin/run_longreads.sh
RUN chmod +x /usr/local/bin/run_longreads.sh
