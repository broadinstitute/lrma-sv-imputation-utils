# Rust helper tools — reference documentation

Three small Rust binaries support the AoU Phase 2 long-reads panel-creation and imputation
pipeline. They are deliberately narrow, streaming, HTSlib-based tools. This directory
documents what each one does, how it is invoked, and the exact math it performs, so that
engineers who did not write them can maintain the pipeline that calls them.

| Tool | Source | One-line purpose |
|---|---|---|
| [`extract_bubble_PLs`](extract-bubble-PLs.md) | `resources/extract-bubble-PLs/` | Re-express an input callset's genotype likelihoods (`PL`) on the reference panel's bi-allelic bubble grid, for GLIMPSE2. |
| [`pop-glimpse2-joint-opt`](pop-glimpse2-joint-opt.md) | `resources/pop-glimpse2/` | Project GLIMPSE2 phased per-path posteriors back onto the constituent atomic bi-allelic variants (`GT:DS:GP`). |
| [`paste-vcfs`](paste-vcfs.md) | `resources/paste-vcfs/` | Column-bind per-sample data from row-aligned cohort shards into one wide BCF. |

## Where they sit in the pipeline

```
 reference panel (bi-allelic bubbles)
        │
        │  input gVCF / joint VCF
        ▼
 [extract_bubble_PLs] ── per-sample PLs on the panel grid ──► GLIMPSE2 (likelihoods → phasing → imputation)
                                                                      │
                                              multi-allelic per-path posteriors (GT:GP)
                                                                      ▼
                                                        [pop-glimpse2-joint-opt] ── per-atomic-variant GT:DS:GP
                                                                      │
                                         (cohort was sharded by sample for the steps above)
                                                                      ▼
                                                            [paste-vcfs] ── single wide multi-sample BCF
```

## Testing and CI

Each tool has an automated, exact-match test suite under
[`resources/tests/`](../resources/tests/). `extract_bubble_PLs` and `paste-vcfs` are
verified against deterministically-generated expected output (integer / structural, so
bit-exact); `pop-glimpse2-joint-opt` uses characterization ("golden") tests because its
output is `f32`. See [`resources/tests/README.md`](../resources/tests/README.md) and the
top-level [`README.md`](../README.md).

## Reference `Cargo.toml` manifests

The three crates already ship their own `Cargo.toml` in the repo. Copies of the
dependency surface each one relies on (reconstructed from the source) are kept under
[`reference/`](reference/) for convenience if a manifest is ever lost or needs to be
recreated. The authoritative manifests are the ones committed next to each crate's `src/`.
