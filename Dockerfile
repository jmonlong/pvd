FROM ubuntu:24.04

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    gcc \
    g++ \
    cmake \
    pkg-config \
    build-essential \
    libhts-dev \
    # libtabixpp-dev \
    # libtabixpp0 \
    tabix \
    apt-transport-https \
    software-properties-common \
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

## bcftools

WORKDIR /build
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
        tar -xjf bcftools-1.20.tar.bz2 && \
        cd bcftools-1.20 && \
        ./configure && make && make install && \
        cd .. && rm -rf bcftools-1.20.tar.bz2

## vcfbub
WORKDIR /bin

RUN wget https://github.com/pangenome/vcfbub/releases/download/v0.1.1/vcfbub && \
    chmod a+x vcfbub

## pvd
WORKDIR /build/pvd

COPY . .

RUN mkdir -p build && cd build && cmake .. && make && mv ../bin/pvd /bin/



WORKDIR /home
