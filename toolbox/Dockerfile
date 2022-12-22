ARG samtoolsversion=1.14
ARG sratoolsversion=3.0.0
ARG reference=NC_045512.2

FROM python:3.7 as final
FROM python:3.7 as build

ARG samtoolsversion
ARG sratoolsversion
ARG reference

RUN apt-get update && apt-get install -qqy python2 wget curl gawk zip zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev procps \
    make python3-pip libxml2-utils bedtools sqlite3 muscle build-essential libssl-dev libcap-dev libexpat1-dev libdbd-sqlite3-perl libboost-dev \
    libboost-program-options-dev libboost-iostreams-dev libboost-timer-dev libboost-chrono-dev libboost-system-dev cpanminus perl-doc vim 

WORKDIR /opt

### samtools
RUN curl -s -L -o htslib-${samtoolsversion}.tar.bz2 https://github.com/samtools/htslib/releases/download/${samtoolsversion}/htslib-${samtoolsversion}.tar.bz2 && \
    tar xjf htslib-${samtoolsversion}.tar.bz2 && rm -f htslib-${samtoolsversion}.tar.bz2 && \
    cd htslib-${samtoolsversion} && \
    ./configure >/dev/null && \
    make --quiet && make --quiet install

RUN curl -s -L -o bcftools-${samtoolsversion}.tar.bz2 https://github.com/samtools/bcftools/releases/download/${samtoolsversion}/bcftools-${samtoolsversion}.tar.bz2 && \
    tar xjf bcftools-${samtoolsversion}.tar.bz2 && \
    cd bcftools-${samtoolsversion} && \
    ./configure --without-curses >/dev/null && \
    make --quiet && make --quiet install

RUN curl -s -L -o samtools-${samtoolsversion}.tar.bz2 https://github.com/samtools/samtools/releases/download/${samtoolsversion}/samtools-${samtoolsversion}.tar.bz2 && \
    tar xjf samtools-${samtoolsversion}.tar.bz2 && rm -f samtools-${samtoolsversion}.tar.bz2 && \
    cd samtools-${samtoolsversion} && \
    ./configure --without-curses >/dev/null && \
    make --quiet && make --quiet install

### mummer
RUN curl -s -L -o mummer-4.0.0rc1.tar.gz https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz && \
    tar xzf mummer-4.0.0rc1.tar.gz && rm mummer-4.0.0rc1.tar.gz && \
    cd mummer-4.0.0rc1 && ./configure >/dev/null && make --quiet && make --quiet install

FROM final

ARG samtoolsversion
ARG sratoolsversion
ARG reference

COPY --from=build /usr/local/bin/ /usr/local/bin
COPY --from=build /usr/local/lib/*mummer* /usr/local/lib/

RUN mkdir -p /pipelines/toolbox/bin /pipelines/toolbox/Scripts/ /usr/local/snpEff /usr/local/trimmomatic

### snpEff
WORKDIR /usr/local/snpEff

RUN curl -s -L -o snpEff_core.zip https://sourceforge.net/projects/snpeff/files/snpEff_v4_2_core.zip/download && \
    unzip -qq snpEff_core.zip && rm -f snpEff_core.zip && mv snpEff 4.2

### Trimmomatic
WORKDIR /usr/local/trimmomatic
RUN curl -s -L -o Trimmomatic-0.39.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && rm -f Trimmomatic-0.39.zip && mv Trimmomatic-0.39 0.39

WORKDIR /opt

RUN apt-get update && apt-get install -qqy python2 wget curl gawk zip procps python3-pip bedtools cpanminus perl-doc vim default-jdk efetch snakemake

### minimap2, The MIT License
RUN curl -s -L -o minimap2-2.17_x64-linux.tar.bz2 https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && \
    tar -jxf minimap2-2.17_x64-linux.tar.bz2 && rm minimap2-2.17_x64-linux.tar.bz2

### SRA-toolkit
RUN curl -s -L -o sratoolkit.${sratoolsversion}-ubuntu64.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${sratoolsversion}/sratoolkit.${sratoolsversion}-ubuntu64.tar.gz && \
    tar xf sratoolkit.${sratoolsversion}-ubuntu64.tar.gz

COPY Dockerfile.user-settings.mkfg /root/.ncbi/user-settings.mkfg
RUN apt-get install uuid-runtime; printf '/LIBS/IMAGE_GUID = "%s"\n' `uuidgen` >> /root/.ncbi/user-settings.mkfg

### Picard
RUN curl -s -L -o picard.jar https://github.com/broadinstitute/picard/releases/download/2.23.9/picard.jar

### GATK
RUN curl -s -L -o gatk-4.2.5.0.zip https://github.com/broadinstitute/gatk/releases/download/4.2.5.0/gatk-4.2.5.0.zip && \
    unzip gatk-4.2.5.0.zip && rm -f gatk-4.2.5.0.zip

RUN python3 -m pip install --upgrade pip==22.2.2
RUN pip3 install Cython protobuf==3.20.0 \
                 jsonpickle pyOpenSSL \
                 biopython nanofilt medaka

RUN curl -s -L -o hisat2-2.2.1-Linux_x86_64.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
    unzip -qq hisat2-2.2.1-Linux_x86_64.zip && rm -f hisat2-2.2.1-Linux_x86_64.zip

ENV PATH=/opt/minimap2-2.17_x64-linux:/opt/sratoolkit.${sratoolsversion}-ubuntu64/bin:/opt/gatk-4.2.5.0:/opt/hisat2-2.2.1:$PATH

RUN PERL_MM_USE_DEFAULT=1; cpan JSON

RUN curl -s -L -o vcftools_0.1.13.tar.gz https://sourceforge.net/projects/vcftools/files/vcftools_0.1.13.tar.gz/download && \
    tar xzf vcftools_0.1.13.tar.gz && rm -f vcftools_0.1.13.tar.gz && cd vcftools_0.1.13 && make

RUN curl -s -L -o mummer2vcf.py https://raw.githubusercontent.com/MatteoSchiavinato/Utilities/master/my-mummer-2-vcf.py && \
    chmod u+x mummer2vcf.py && mv mummer2vcf.py /pipelines/toolbox/Scripts/

RUN curl -s -L -o vcf_validator_linux https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.4/vcf_validator_linux && \
    chmod u+x vcf_validator_linux

RUN ln -s /bin/bash /usr/bin/bash
RUN ln -s /usr/bin/python2 /usr/bin/python

### Bio:DB::SAM
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y bioperl libncurses5-dev jq

RUN curl -s -L -o Bio-SamTools-1.43.tar.gz https://cpan.metacpan.org/authors/id/L/LD/LDS/Bio-SamTools-1.43.tar.gz && \
    tar xzf Bio-SamTools-1.43.tar.gz && rm -f Bio-SamTools-1.43.tar.gz && \
    cd Bio-SamTools-1.43 && chmod 600 INSTALL.pl && sed -i 's/sudo//g' INSTALL.pl && (echo y|perl INSTALL.pl)

WORKDIR /pipelines/toolbox/

COPY . /pipelines/toolbox/

RUN gatk CreateSequenceDictionary -R /pipelines/toolbox/static/reference/${reference}.fa && \
    samtools faidx /pipelines/toolbox/static/reference/${reference}.fa && \
    gatk IndexFeatureFile --input /pipelines/toolbox/static/reference/${reference}.known_sites.vcf && \
    hisat2-build /pipelines/toolbox/static/reference/${reference}.fa /pipelines/toolbox/static/reference/${reference}.fa

RUN medaka tools download_models --model r941_prom_high_g360 r941_min_high_g360
RUN sed -i 's/python/python3/g' /usr/local/bin/medaka_consensus

ENV HISAT2_HOME=/opt/hisat2-2.2.1 \
    PERL5LIB=/opt/vcftools_0.1.13/perl/ \
    LD_LIBRARY_PATH=/usr/local/lib \
    PATH=/opt:/opt/vcftools_0.1.13/bin:/usr/local/snpEff/4.2/scripts:/pipelines/toolbox/Scripts:$PATH

WORKDIR /pipelines/toolbox/workflow
