# syntax=docker/dockerfile:1
#
# Build environment and compiled binaries for the code in this repository.
# Tested environment of the paper: Ubuntu 24.04, g++ 13.3.0, Armadillo 14.4.0,
# GSL 2.8, Boost 1.83.0, pthread OpenBLAS (BLAS + LAPACK).
#
# Build:  docker build -t rep .
# Run:    see README.md ("Docker" section)

# TODO (optional, for stricter reproducibility): pin the base image by digest, e.g.
#   FROM ubuntu:24.04@sha256:<digest>
# Get the digest with: docker pull ubuntu:24.04 && docker inspect --format='{{index .RepoDigests 0}}' ubuntu:24.04
FROM ubuntu:24.04@sha256:73b999393d07f2bf65a40373778a44dc9a6d54ad7718d29247d8ff9413c4d6a8

ARG GSL_VERSION=2.8
ARG ARMA_VERSION=14.4.0
ARG GSL_SHA256="6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190"
ARG ARMA_SHA256="023242fd59071d98c75fb015fd3293c921132dc39bf46d221d4b059aae8d79f4"

ENV DEBIAN_FRONTEND=noninteractive

# ---------------------------------------------------------------------------
# 1. Toolchain, Boost and OpenBLAS/LAPACK from the Ubuntu 24.04 repositories
#    (Ubuntu 24.04 ships g++ 13.3.0 and Boost 1.83.0)
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        xz-utils \
        libboost-dev \
        libboost-program-options-dev \
        libopenblas-pthread-dev \
        liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. GSL (built from source)
# ---------------------------------------------------------------------------
WORKDIR /tmp/build
RUN curl -fsSL -o gsl.tar.gz "https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz" \
    && if [ -n "${GSL_SHA256}" ]; then echo "${GSL_SHA256}  gsl.tar.gz" | sha256sum -c -; fi \
    && tar xzf gsl.tar.gz \
    && cd "gsl-${GSL_VERSION}" \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd .. && rm -rf gsl.tar.gz "gsl-${GSL_VERSION}"

# ---------------------------------------------------------------------------
# 3. Armadillo (built from source)
# ---------------------------------------------------------------------------
RUN curl -fsSL -o arma.tar.xz "https://sourceforge.net/projects/arma/files/armadillo-${ARMA_VERSION}.tar.xz/download" \
    && if [ -n "${ARMA_SHA256}" ]; then echo "${ARMA_SHA256}  arma.tar.xz" | sha256sum -c -; fi \
    && tar xf arma.tar.xz \
    && cd "armadillo-${ARMA_VERSION}" \
    && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build build -j"$(nproc)" \
    && cmake --install build \
    && cd .. && rm -rf arma.tar.xz "armadillo-${ARMA_VERSION}"

RUN ldconfig

# ---------------------------------------------------------------------------
# 4. Print the versions that end up in the image (visible in the build log)
# ---------------------------------------------------------------------------
RUN g++ --version | head -n 1 \
    && dpkg-query -W -f='Boost (apt): ${Version}\n' libboost-dev \
    && ldconfig -p | grep -E 'libarmadillo|libgsl' \
    && (update-alternatives --display libblas.so.3 | head -n 2 || true) \
    && (update-alternatives --display liblapack.so.3 | head -n 2 || true)

# ---------------------------------------------------------------------------
# 5. Repository code and compilation (same commands as in the README)
# ---------------------------------------------------------------------------
WORKDIR /rep
COPY . /rep

RUN g++ -O3 -o RNN_only.out RNN_only.cpp -lgsl -lgslcblas -lm -larmadillo -lboost_program_options -fopenmp \
    && g++ -O3 -o "RNN+FFL.out" "RNN+FFL.cpp" -lgsl -lgslcblas -lm -larmadillo -lboost_program_options -fopenmp \
    && mkdir -p simulation_results

# Default number of OpenMP threads (as in the paper); override with: docker run -e OMP_NUM_THREADS=N ...
ENV OMP_NUM_THREADS=6

CMD ["bash"]
