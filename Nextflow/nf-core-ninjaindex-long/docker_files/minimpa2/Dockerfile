# Base Image
FROM nfcore/base
LABEL authors="Xiandong Meng" \
      description="Docker image containing all requirements for nf-core/ninjaindex pipeline" \
      software.name="2" \
      software.description="A long reads mapping tool" \
      software.website="https://github.com/lh3/minimap2" \
      software.documentation="https://github.com/lh3/minimap2" \
      container.category="tool"

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install minimap2
RUN conda install -c bioconda -y minimap2 samtools bbmap
RUN conda install -c conda-forge -y awscli

ADD run_minimap2.sh /usr/local/bin/run_minimap2.sh
RUN chmod +x /usr/local/bin/run_minimap2.sh
