FROM ubuntu:26.04

LABEL maintainer="Pipeline Maintainer"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# ----------------------------------------------------------------------
# 1. System & Bioinformatics Dependencies (APT)
# ----------------------------------------------------------------------
# Enable 32-bit architecture for the original 32-bit Muscle
RUN dpkg --add-architecture i386 \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    # 32-bit compatibility libraries
    libc6:i386 \
    libstdc++6:i386 \
    # Build & System Tools
    apt-utils \
    automake \
    autoconf \
    autotools-dev \
    build-essential \
    bzip2 \
    ca-certificates \
    curl \
    git \
    wget \
    # Libraries
    expat \
    libexpat1-dev \
    libgd-dev \
    libglib2.0-dev \
    libidn12 \
    # Perl & Modules
    perl \
    bioperl \
    libgd-perl \
    liblist-moreutils-perl \
    libwww-perl \
    libxml-dom-perl \
    libxml-dom-xpath-perl \
    # Bioinformatics Packages available via APT
    aragorn \
    emboss \
    exonerate \
    hmmer \
    infernal \
    mafft \
    ncbi-blast+ \
    # Python Environment
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# ----------------------------------------------------------------------
# 2. Binary Tools (Original 32-bit Muscle, tbl2asn, Erpin)
# ----------------------------------------------------------------------
WORKDIR /tmp/tools

# Original 32-bit Muscle (does not segfault on modern Linux)
RUN wget -q http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz \
    && tar -xzf muscle3.8.31_i86linux32.tar.gz \
    && mv muscle3.8.31_i86linux32 /usr/local/bin/muscle \
    && chmod +x /usr/local/bin/muscle \
    && rm -f muscle3.8.31_i86linux32.tar.gz

# tbl2asn
RUN wget -q https://anaconda.org/bioconda/tbl2asn/25.7/download/linux-64/tbl2asn-25.7-0.tar.bz2 \
    && tar -xjf tbl2asn-25.7-0.tar.bz2 \
    && chmod 755 bin/tbl2asn \
    && mv bin/tbl2asn /usr/local/bin/tbl2asn \
    && rm -rf /tmp/tools/* \
    && ln -s /usr/lib/x86_64-linux-gnu/libidn.so.12 /usr/lib/x86_64-linux-gnu/libidn.so.11

# Erpin
   RUN wget -q http://rssf.i2bc.paris-saclay.fr/download/Erpin/erpin5.5.4.serv.tar.gz \
    && tar -xzf erpin5.5.4.serv.tar.gz \
    && cp erpin5.5.4.serv/bin/erpin /usr/local/bin/ \
    && chmod +x /usr/local/bin/erpin \
    && rm -rf /tmp/tools/*

# ----------------------------------------------------------------------
# 3. BFL-lab Tools & Repositories
# ----------------------------------------------------------------------
# Trust all git directories globally to avoid "dubious ownership" fatal errors
RUN git config --global --add safe.directory "*"

# PirObject & PirModels
RUN git clone --depth 1 https://github.com/prioux/PirObject.git /tmp/PirObject \
    && cp /tmp/PirObject/lib/PirObject.pm /etc/perl/ \
    && rm -rf /tmp/PirObject \
    && git clone --depth 1 https://github.com/BFL-lab/PirModels.git /PirModels

# flip
RUN git clone --depth 1 https://github.com/BFL-lab/flip.git /tmp/flip \
    && gcc -o /usr/local/bin/flip /tmp/flip/src/flip.c \
    && chmod +x /usr/local/bin/flip \
    && rm -rf /tmp/flip

# umac, HMMsearchWC, CMsearchW
RUN git clone --depth 1 https://github.com/BFL-lab/umac.git /tmp/umac \
    && cp /tmp/umac/umac /usr/local/bin/ \
    && git clone --depth 1 https://github.com/BFL-lab/HMMsearchWC.git /tmp/HMMsearchWC \
    && cp /tmp/HMMsearchWC/HMMsearchCombiner /tmp/HMMsearchWC/HMMsearchWrapper /usr/local/bin/ \
    && git clone --depth 1 https://github.com/BFL-lab/CMsearchW.git /tmp/CMsearchW \
    && cp /tmp/CMsearchW/CMsearchWrapper /usr/local/bin/ \
    && rm -rf /tmp/umac /tmp/HMMsearchWC /tmp/CMsearchW

# RNAfinder
RUN git clone --depth 1 https://github.com/BFL-lab/RNAfinder.git /tmp/RNAfinder \
    && cp /tmp/RNAfinder/RNAfinder /usr/local/bin/ \
    && cp /tmp/RNAfinder/DOT_RNAfinder.cfg /.RNAfinder.cfg \
    && rm -rf /tmp/RNAfinder

# grab-fasta
RUN git clone --depth 1 https://github.com/BFL-lab/grab-fasta.git /tmp/grab-fasta \
    && cp /tmp/grab-fasta/grab-fasta /tmp/grab-fasta/grab-seq /usr/local/bin/ \
    && rm -rf /tmp/grab-fasta

# mf2sqn
RUN git clone --depth 1 https://github.com/BFL-lab/mf2sqn.git /mf2sqn \
    && cp /mf2sqn/mf2sqn /usr/local/bin/ \
    && mkdir -p /usr/share/perl5 \
    && cp /mf2sqn/qualifs.pl /usr/share/perl5/

# MFannot & Data
RUN git clone https://github.com/BFL-lab/mfannot.git /mfannot \
    && ln -s /mfannot/mfannot /usr/local/bin/mfannot \
    && cp -r /mfannot/examples /examples \
    && git clone --depth 1 https://github.com/BFL-lab/MFannot_data.git /MFannot_data

# Ensure all scripts/binaries in /usr/local/bin have execute permissions
RUN chmod +x /usr/local/bin/* /mfannot/mfannot

# BLAST matrices
RUN mkdir -p /BLASTMAT \
    && ( cd /BLASTMAT && wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/* || \
         wget -q -e robots=off -r -np -nd https://ftp.ncbi.nlm.nih.gov/blast/matrices/ -P /BLASTMAT/ ) \
    && test -f /BLASTMAT/PAM70

# ----------------------------------------------------------------------
# 4. Environment Variables
# ----------------------------------------------------------------------
ENV RNAFINDER_CFG_PATH="/" \
    MF2SQN_LIB="/mf2sqn/lib/" \
    MFANNOT_LIB_PATH="/MFannot_data/protein_collections/" \
    MFANNOT_EXT_CFG_PATH="/MFannot_data/config" \
    MFANNOT_MOD_PATH="/MFannot_data/models/" \
    BLASTMAT="/BLASTMAT/" \
    EGC="/MFannot_data/EGC/" \
    ERPIN_MOD_PATH="/MFannot_data/models/Erpin_models/" \
    PIR_DATAMODEL_PATH="/PirModels" \
    PERL5LIB="/mf2sqn/lib" \
    PATH="/opt/venv/bin:/mfannot:${PATH}"

# ----------------------------------------------------------------------
# 5. Pipeline Python Setup & Application Code
# ----------------------------------------------------------------------
RUN python3 -m venv /opt/venv

WORKDIR /opt/pipeline

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY src/ /opt/pipeline/src/
WORKDIR /opt/pipeline/src

ENTRYPOINT ["python3", "run.py"]