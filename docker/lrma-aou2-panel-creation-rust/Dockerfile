# ==========================================
# Stage 1: Optimized bcftools Builder
# ==========================================
FROM ubuntu:24.04 AS builder

# Install dependencies needed for an optimal bcftools build
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdeflate-dev \
    libgsl-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src

# Download, compile, and stage the latest BCFtools with full optimizations
ENV BCFTOOLS_VERSION=1.23.1
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -vxjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure \
    --prefix=/usr/local \
    --enable-libcurl \
    --with-libdeflate \
    --enable-gsl && \
    make -j$(nproc) && \
    make install DESTDIR=/tmp/bcftools_stage

# ==========================================
# Stage 2: Development & Runtime Environment
# ==========================================
FROM ubuntu:24.04

# 1. Install build tools + dev headers for Rust/HTSlib compilation
# 2. Add libdeflate-dev & libgsl-dev for runtime execution
# 3. Add perl (for bcftools scripts)
# 4. Install Google Cloud CLI
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    clang \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdeflate-dev \
    libgsl-dev \
    perl \
    apt-transport-https \
    ca-certificates \
    gnupg \
    curl \
    python3 \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    google-cloud-cli \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Install Rust compiler and cargo
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Copy the optimized bcftools (including plugins) from the staging directory
COPY --from=builder /tmp/bcftools_stage/usr/local /usr/local/

# Explicitly declare the plugin directory for bcftools
ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"
ENV HOME=/root

# Set the default entrypoint
CMD ["bash"]
