# lrma-aou2 Rust tools — build, test & package

A lightweight, reproducible CI/CD setup for the three Rust helper tools used in the
`broadinstitute/lrma-aou2-panel-creation` pipeline (branch `sl_aou2_v1`):

| Tool | Binary | One-line purpose |
|---|---|---|
| `extract-bubble-PLs` | `extract-bubble-PLs` | Project PLs from a gVCF/joint VCF onto a panel's bi-allelic "bubble" sites (or a subset thereof) for GLIMPSE2. |
| `pop-glimpse2` | `pop-glimpse2` | Project GLIMPSE2 phased per-path posteriors back onto the constituent atomic bi-allelic variants (yielding `GT:DS:GP`). |
| `paste-vcfs` | `paste-vcfs` | Essentially a lean version of bcftools merge; horizontally concatenate sample columns across VCF/BCFs that share identical sites. |

Per-tool design, every parameter, and the exact math are documented under [`docs/`](docs/).

## Layout

```
rust-toolchain.toml            pinned toolchain (rustc >= 1.86; see below)
resources/<tool>/              each tool is an INDEPENDENT crate:
    Cargo.toml                   its manifest (upstream, verbatim)
    Cargo.lock                   its committed, known-good lockfile
    src/                         its sources
resources/tests/               exact-match integration suite (python3 + bcftools only)
docker/lrma-aou2-panel-creation-rust/Dockerfile   verbatim upstream base image (bcftools + rust)
docker/app/Dockerfile          multi-stage build of the three release binaries
docs/                          maintainer documentation, one file per tool
.github/workflows/             ci.yml (native build+test) and docker-image.yml (test inside image)
Makefile                       thin wrappers around the exact commands CI runs
```

The three tools are **separate crates, not a Cargo workspace** — each has its own
`Cargo.lock`, exactly as they were built upstream. This matters because their dependency
graphs differ (only two use `rust-htslib`) and they even use different editions
(`pop-glimpse2` is edition 2024; the others are 2021), which a single shared lock could not
represent faithfully.

## Reproducibility: committed lockfiles + pinned toolchain

Each crate's `resources/<tool>/Cargo.lock` is committed and every build uses `--locked`, so
dependency versions are frozen to the exact set that was validated upstream (e.g.
`rust-htslib 0.44.1`, `hts-sys 2.2.0`, `clap 4.6.1`). Nothing needs generating; if you ever
deliberately update dependencies, run `make lock` (per-crate `cargo generate-lockfile`) and
commit the changed locks.

The toolchain is pinned to **Rust 1.90** in `rust-toolchain.toml` — the version the original
builds used. The hard floor is 1.86 regardless: `pop-glimpse2` is `edition = "2024"` (needs
rustc ≥ 1.85), and the `extract`/`paste` locks pin `icu 2.2.0` (pulled in via
`rust-htslib` → `url` → `idna`'s default ICU4X back end), which needs rustc ≥ 1.86. `rustup`
fetches the pinned version automatically. Any newer stable also works; building on an older
toolchain will fail while parsing/resolving those dependencies, which is a toolchain-floor
issue rather than a problem with this repo.

## Building and testing locally

```bash
make build      # build each crate --release --locked (per-crate target dirs)
make test       # build, collect binaries into ./dist, run the exact-match suite
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
  the cargo cache, then builds each crate `--release --locked` against its committed lock,
  collects the three binaries into `./dist`, and runs the exact-match suite. There is **no
  formatting or lint gate** — `cargo build` + the suite are the only checks, so style/clippy
  never blocks a merge.
* **`docker-image.yml`** (on demand / release branches): builds the base image and the app
  image with plain `docker build` (so the app's `FROM` sees the locally-built base), then
  runs the suite inside the final image.

## Integrating into the upstream repo

This scaffold **vendors** the three tool sources under `resources/<tool>/src` (transcribed
from `sl_aou2_v1`) so the suite is self-contained and CI can compile them. When merging
into the real repository:

* Each tool stays an independent crate with its own `Cargo.toml` + `Cargo.lock`. If the
  upstream tree already has these crates, delete the vendored `src/` copies here and point
  the build/CI/Docker steps at the existing crate directories (they only need the
  `resources/<tool>/` paths and the committed locks).
* Keep the per-crate `Cargo.lock` files under version control — they are the reproducibility
  contract that this setup relies on.
* The directory paths here mirror upstream (`resources/<tool>/…`,
  `docker/lrma-aou2-panel-creation-rust/…`), so `docs/`, `resources/tests/`, `docker/app/`,
  the workflows, and `rust-toolchain.toml` can be dropped in with minimal reshuffling.
