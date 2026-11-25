# oligoscan + bowtie1 + bowtie2 + RIsearch2 in one image.
#   bowtie1    -> exhaustive oracle for validating oligoscan at <= 3 MM
#   bowtie2    -> comparison runs against the current platform
#   RIsearch2  -> energy-based seed-and-extend comparator (Alkan 2017, NAR)
#   oligoscan itself needs NO index; it reads the FASTA directly
FROM ubuntu:24.04

# -march=native targets the CPU that builds the image. Fine for a local image;
# set MARCH=x86-64-v3 (AVX2) or x86-64-v2 for an image that runs elsewhere.
ARG MARCH=native

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        bowtie bowtie2 \
        build-essential \
        python3 python3-numpy python3-pip \
        pigz ca-certificates wget curl bzip2 time \
        zlib1g libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# ViennaRNA (Turner 2004 duplexfold) + pyfaidx for the case-study pipeline.
# Same library call the siRNA/ASO designer app uses, so energies agree; kept
# here so the case study is reproducible from this repository alone.
RUN pip install --no-cache-dir --break-system-packages \
        "ViennaRNA>=2.6" "pyfaidx>=0.8" "numpy>=1.26" \
    && python3 -c "import RNA, pyfaidx; print('ViennaRNA', RNA.__version__)"

WORKDIR /opt/oligoscan
COPY src/     src/
COPY tests/   tests/
COPY scripts/ scripts/
RUN chmod +x scripts/*.sh && ln -s /opt/oligoscan/scripts/*.sh /usr/local/bin/

# Build, then refuse to produce an image whose numerics or selection layer fail.
RUN gcc -O3 -march=${MARCH} -fopenmp -Wall -o /usr/local/bin/oligoscan src/oligoscan.c \
    && oligoscan --selftest \
    && python3 tests/test_oligoscan_select.py /usr/local/bin/oligoscan \
    && python3 tests/test_junction_edge.py --bin /usr/local/bin/oligoscan \
    && python3 tests/test_bulge_leading_mismatch.py --bin /usr/local/bin/oligoscan \
    && python3 tests/test_exhaustiveness.py --bin /usr/local/bin/oligoscan --quick

# RIsearch2 v2.1 -- the PUBLISHED C tool (Alkan et al., NAR 2017, 45(8):e60),
# GPL3, taken from the authors' own conda channel and pinned by md5.
# NOT github.com/RTH-tools/risearch, which now hosts an unrelated Rust rewrite
# (3.0.0-alpha, BUSL-1.1, citation in preparation).
ARG RISEARCH_PKG=risearch-2.1-h9bf148f_0.tar.bz2
ARG RISEARCH_MD5=82b021a622ef6dc9ea187242c0154e9c
RUN cd /tmp \
    && curl -fsSL -o ris.tar.bz2 "https://conda.anaconda.org/rthtools/linux-64/${RISEARCH_PKG}" \
    && echo "${RISEARCH_MD5}  ris.tar.bz2" | md5sum -c - \
    && mkdir -p ris && tar -xjf ris.tar.bz2 -C ris \
    && install -m 0755 ris/bin/RIsearch2 /usr/local/bin/RIsearch2 \
    && install -m 0755 ris/bin/RIsearch1 /usr/local/bin/RIsearch1 \
    && rm -rf /tmp/ris /tmp/ris.tar.bz2 \
    && RIsearch2 -h 2>&1 | grep -q 'RIsearch v2.1'

RUN echo "bowtie1: $(bowtie --version | head -1)" > /opt/oligoscan/VERSIONS \
    && echo "bowtie2: $(bowtie2 --version | head -1)" >> /opt/oligoscan/VERSIONS \
    && echo "RIsearch2: $(RIsearch2 -h 2>&1 | grep -o 'RIsearch v[0-9.]*' | head -1) (conda rthtools, md5 ${RISEARCH_MD5})" >> /opt/oligoscan/VERSIONS \
    && echo "ViennaRNA: $(python3 -c 'import RNA; print(RNA.__version__)')" >> /opt/oligoscan/VERSIONS

WORKDIR /data
CMD ["bash"]
