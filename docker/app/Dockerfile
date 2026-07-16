# syntax=docker/dockerfile:1.7
# Final image containing the three compiled Rust tools.
#
# Built FROM the repository's existing base image
# (docker/lrma-aou2-panel-creation-rust/Dockerfile), which already provides Rust
# + cargo, the HTSlib system dependencies that rust-htslib needs to build its
# vendored HTSlib, and a from-source bcftools at /usr/local/bin/bcftools.
#
# Build it, pointing BASE_IMAGE at a prebuilt tag to skip recompiling bcftools:
#   docker build -f docker/app/Dockerfile \
#       --build-arg BASE_IMAGE=lrma-aou2-panel-creation-rust:latest \
#       -t lrma-sv-imputation-utils:latest .
#
# The three tools are INDEPENDENT crates, each with its own committed Cargo.lock,
# and are built separately with `--locked` (reproducing the original builds; see
# rust-toolchain.toml for the required rustc >= 1.86). The toolchain file pins
# the exact rustc, which rustup fetches on first cargo invocation.
#
# The final stage copies ONLY the stripped binaries. No test code or fixtures
# are ever added to any layer (see .dockerignore); the suite is bind-mounted at
# run time instead.

ARG BASE_IMAGE=lrma-aou2-panel-creation-rust:latest

# ---------------------------------------------------------------------------
# Stage 1: build each crate on the base toolchain image, against its own lock.
# ---------------------------------------------------------------------------
FROM ${BASE_IMAGE} AS build
WORKDIR /src

# Install the pinned toolchain first (cheap layer, cached across source churn).
COPY rust-toolchain.toml ./
RUN rustup show >/dev/null 2>&1 || true

# Copy each crate (manifest + lock + sources) and build it with --locked.
COPY resources/extract-bubble-PLs resources/extract-bubble-PLs
COPY resources/pop-glimpse2       resources/pop-glimpse2
COPY resources/paste-vcfs         resources/paste-vcfs

RUN --mount=type=cache,target=/root/.cargo/registry \
    set -eux; \
    for c in extract-bubble-PLs pop-glimpse2 paste-vcfs; do \
        cargo build --release --locked --manifest-path resources/$c/Cargo.toml; \
    done; \
    mkdir -p /out; \
    cp resources/extract-bubble-PLs/target/release/extract-bubble-PLs /out/; \
    cp resources/pop-glimpse2/target/release/pop-glimpse2             /out/; \
    cp resources/paste-vcfs/target/release/paste-vcfs                 /out/; \
    strip /out/* || true

# ---------------------------------------------------------------------------
# Stage 2: lean runtime. Same base (for bcftools + shared libs); only the three
# binaries are added on top -- no Rust sources, no cargo caches, no tests.
# ---------------------------------------------------------------------------
FROM ${BASE_IMAGE}

COPY --from=build /out/extract-bubble-PLs      /usr/local/bin/extract-bubble-PLs
COPY --from=build /out/pop-glimpse2            /usr/local/bin/pop-glimpse2
COPY --from=build /out/paste-vcfs              /usr/local/bin/paste-vcfs

# Sanity check that a binary is on PATH and dynamically resolvable.
RUN paste-vcfs --help >/dev/null 2>&1 || true

CMD ["bash"]
