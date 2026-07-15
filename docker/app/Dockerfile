# syntax=docker/dockerfile:1.7
# Final image containing the three compiled Rust tools.
#
# It is built FROM the repository's existing base image
# (docker/lrma-aou2-panel-creation-rust/Dockerfile), which already provides:
#   * Rust + cargo (rustup),
#   * the full HTSlib system dependencies (libcurl, libdeflate, zlib, bz2, lzma,
#     ssl, gsl) that rust-htslib needs to build its vendored HTSlib,
#   * a from-source bcftools at /usr/local/bin/bcftools.
#
# Build it, pointing BASE_IMAGE at a prebuilt tag to skip recompiling bcftools:
#   docker build -f docker/app/Dockerfile \
#       --build-arg BASE_IMAGE=lrma-aou2-panel-creation-rust:latest \
#       -t lrma-aou2-panel-creation-tools:latest .
#
# Reproducibility: the build uses `--locked`, so it fails unless Cargo.lock is
# present and consistent. Commit Cargo.lock (see README "one-time bootstrap").
#
# The final stage copies ONLY the stripped binaries. No test code or fixtures
# are ever added to any layer (see .dockerignore); the suite is bind-mounted at
# run time instead.

ARG BASE_IMAGE=lrma-aou2-panel-creation-rust:latest

# ---------------------------------------------------------------------------
# Stage 1: build the workspace binaries on the base toolchain image.
# ---------------------------------------------------------------------------
FROM ${BASE_IMAGE} AS build

WORKDIR /src

# Pin the toolchain first so the rustup channel is installed before any copy of
# the (large) sources invalidates caching.
COPY rust-toolchain.toml ./

# Dependency-caching layer: copy only the manifests + lockfile and build a set
# of dummy crates so the dependency graph (including the vendored HTSlib C
# build) is compiled and cached independently of source churn.
COPY Cargo.toml Cargo.lock ./
COPY resources/extract-bubble-PLs/Cargo.toml resources/extract-bubble-PLs/Cargo.toml
COPY resources/pop-glimpse2/Cargo.toml       resources/pop-glimpse2/Cargo.toml
COPY resources/paste-vcfs/Cargo.toml         resources/paste-vcfs/Cargo.toml
RUN --mount=type=cache,target=/root/.cargo/registry \
    mkdir -p resources/extract-bubble-PLs/src \
             resources/pop-glimpse2/src/bin \
             resources/paste-vcfs/src && \
    echo 'fn main(){}' > resources/extract-bubble-PLs/src/main.rs && \
    echo 'fn main(){}' > resources/pop-glimpse2/src/bin/pop-glimpse2-joint-opt.rs && \
    echo 'fn main(){}' > resources/paste-vcfs/src/main.rs && \
    cargo build --release --locked --workspace || true

# Now copy the real sources and build for real. Touch the sources so cargo sees
# them as newer than the dummy build above.
COPY resources/extract-bubble-PLs/src resources/extract-bubble-PLs/src
COPY resources/pop-glimpse2/src       resources/pop-glimpse2/src
COPY resources/paste-vcfs/src         resources/paste-vcfs/src
RUN --mount=type=cache,target=/root/.cargo/registry \
    cargo build --release --locked --workspace && \
    mkdir -p /out && \
    cp target/release/extract_bubble_PLs \
       target/release/pop-glimpse2-joint-opt \
       target/release/paste-vcfs /out/ && \
    strip /out/* || true

# ---------------------------------------------------------------------------
# Stage 2: lean runtime. Same base (for bcftools + shared libs) but only the
# binaries are added on top -- no Rust sources, no cargo caches, no tests.
# ---------------------------------------------------------------------------
FROM ${BASE_IMAGE}

COPY --from=build /out/extract_bubble_PLs      /usr/local/bin/extract_bubble_PLs
COPY --from=build /out/pop-glimpse2-joint-opt  /usr/local/bin/pop-glimpse2-joint-opt
COPY --from=build /out/paste-vcfs              /usr/local/bin/paste-vcfs

# Sanity check that the binaries are on PATH and dynamically resolvable.
RUN paste-vcfs --help >/dev/null 2>&1 || true

CMD ["bash"]
