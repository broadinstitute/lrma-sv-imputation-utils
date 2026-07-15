# lrma-aou2 Rust tools — build, test & package

A lightweight, reproducible CI/CD setup for the three Rust helper tools used in the
`broadinstitute/lrma-aou2-panel-creation` pipeline (branch `sl_aou2_v1`):

| Tool | Binary | One-line purpose |
|---|---|---|
| `extract-bubble-PLs` | `extract_bubble_PLs` | Project PLs from a gVCF/joint VCF onto a panel's bi-allelic "bubble" sites for GLIMPSE2. |
| `pop-glimpse2` | `pop-glimpse2-joint-opt` | Project phased GLIMPSE2 joint posteriors from multi-allelic paths back onto atomic bi-allelic variants. |
| `paste-vcfs` | `paste-vcfs` | Horizontally concatenate sample columns across VCF/BCFs that share identical sites. |

Per-tool design, every parameter, and the exact math are documented under [`docs/`](docs/).

## Layout

```
Cargo.toml                     workspace (one Cargo.lock for all three crates)
rust-toolchain.toml            pinned toolchain (rustc + rustfmt + clippy)
resources/<tool>/              each tool's crate (Cargo.toml + src) — mirrors upstream paths
resources/tests/               exact-match integration suite (python3 + bcftools only)
docker/lrma-aou2-panel-creation-rust/Dockerfile   verbatim upstream base image (bcftools + rust)
docker/app/Dockerfile          multi-stage build of the three release binaries
docs/                          maintainer documentation, one file per tool
.github/workflows/             ci.yml (native build+test) and docker-image.yml (test inside image)
Makefile                       thin wrappers around the exact commands CI runs
```

## One-time bootstrap: create and commit `Cargo.lock`

Reproducibility hinges on a committed `Cargo.lock`, and every build here uses `--locked`.
The lock is **not** generated in this scaffold because valid dependency checksums require
network access to crates.io. Run this once, then commit the result:

```bash
cargo generate-lockfile      # or: make lock
git add Cargo.lock
git commit -m "Pin dependency versions"
```

`rust-htslib` is pinned to `0.44` in the crate manifests. If that minor turns out to
predate an API these sources use, bump the `rust-htslib` line in the three
`resources/*/Cargo.toml`, re-run `cargo generate-lockfile`, and commit again. This is the
only manual step; after it, CI and Docker builds are fully deterministic.

## Building and testing locally

```bash
make build      # cargo build --workspace --locked  (debug)
make test       # build, then run the exact-match suite against target/debug
make fmt clippy # formatting + lints, exactly as CI enforces them
```

Or drive the suite directly against any binaries:

```bash
resources/tests/run_all.sh --bin-dir target/release
# or point individual tools somewhere else:
PASTE_VCFS_BIN=/usr/local/bin/paste-vcfs resources/tests/run_all.sh
```

The suite needs only `python3` (standard library) and `bcftools` — the same two tools the
runtime image ships — so it runs unchanged inside the final container.

### How the tests stay exact

* **extract-bubble-PLs** — integer/string logic only, so `generate.py` is a bit-exact
  oracle (it even double-checks itself against a hand-derived expectation). Fixtures are
  turned into indexed BCF/VCF.GZ via `bcftools`, and `bcftools query` output is compared
  byte-for-byte.
* **paste-vcfs** — no arithmetic; the expected output is *derived from the inputs* by
  re-querying each with the same `bcftools` format and horizontally joining, so any float
  rendering cancels out on both sides of the diff.
* **pop-glimpse2** — `f32` arithmetic, so `generate.py` is a faithful `f32`
  re-implementation and picks fixtures by **rejection sampling**: any draw landing near a
  rounding or threshold boundary (GT at 0.5, DS/GP on the 0.001/permille grid, or the
  top-k score cut) is rejected, guaranteeing the Rust `f32` output is byte-identical.

Every generator is deterministic (fixed seeds), so expected outputs are stable across runs.

## Docker

The final image is built **FROM the existing base image**, which already provides Rust and
all HTSlib system dependencies plus a from-source `bcftools`:

```bash
# Build the base once (slow: compiles bcftools from source) …
docker build -f docker/lrma-aou2-panel-creation-rust/Dockerfile \
    -t lrma-aou2-panel-creation-rust:latest .

# … then build the tools image on top (fast), and test the shipped binaries:
make docker      BASE_IMAGE=lrma-aou2-panel-creation-rust:latest
make docker-test BASE_IMAGE=lrma-aou2-panel-creation-rust:latest
```

`make docker-test` runs the exact-match suite **inside** the final image against
`/usr/local/bin`, with `resources/tests/` bind-mounted read-only. Tests and fixtures are
never baked into any image layer (enforced by `.dockerignore`), so the shipped image
contains only the three stripped release binaries on top of the base.

If you publish the base image to a registry, point `BASE_IMAGE` at it (and, in
`docker-image.yml`, replace the base-build step with a `docker pull`) to skip the bcftools
compile and keep CI well within budget.

## CI

* **`ci.yml`** (every PR/push): apt-installs the HTSlib dev headers + `bcftools`, restores
  the cargo cache, checks `fmt`/`clippy`, builds `--locked` (debug), and runs the suite.
  Debug builds keep it under a ~10-minute budget; the tools are deterministic, so
  optimisation level does not affect output.
* **`docker-image.yml`** (on demand / release branches): builds the base + app images with
  layer caching and runs the suite inside the final image.

## Integrating into the upstream repo

This scaffold **vendors** the three tool sources under `resources/<tool>/src` (transcribed
from `sl_aou2_v1`) so the suite is self-contained and CI can compile them. When merging
into the real repository:

* Keep a single copy of each tool's `src/`. If the upstream tree already has these
  sources, delete the vendored copies here and let the workspace `members` point at the
  existing crate directories.
* If the repo root already has a `Cargo.toml`, merge the `members` list from this one into
  it rather than overwriting.
* The directory paths here mirror upstream (`resources/<tool>/…`,
  `docker/lrma-aou2-panel-creation-rust/…`), so `docs/`, `resources/tests/`, `docker/app/`,
  the workflows, and the workspace files can be dropped in with minimal reshuffling.
