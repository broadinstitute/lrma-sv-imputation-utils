# `paste-vcfs`

Source: `resources/paste-vcfs/src/main.rs`

## 1. High-level function

`paste-vcfs` is the VCF analogue of the Unix `paste` command: it takes several
BCF/`VCF.GZ` files that describe **exactly the same sites in exactly the same order**, each
holding a *different set of samples*, and glues their sample columns together side-by-side
into one multi-sample BCF. Selected `INFO` and `FORMAT` fields are carried through.

In one sentence: **it column-binds per-sample data from cohort shards that share an
identical site list, producing one wide BCF.**

This is used to reassemble a cohort that was split by sample for parallel processing (e.g.
per-shard GLIMPSE2 runs) back into a single callset, without the cost of a coordinate merge.

## 2. Design

Correctness rests on the caller's guarantee that all inputs are **row-aligned**: the same
number of records, in the same order, with identical `CHROM/POS/ID/REF/ALT`. The tool reads
one record from every input in lockstep and **asserts** that alignment on every row; any
divergence is a `FATAL` panic rather than a silent mismerge.

Two read paths are unified behind a `Source` enum:

* **Whole-file mode** (no `--region`) uses the safe high-level `bcf::Reader`.
* **Region mode** (`--region`) uses a small `RegionReader` that drives HTSlib's CSI index
  iterator directly through the raw FFI (`hts_idx_load` + `hts_itr_query` +
  `hts_itr_next`). The contig id is resolved once via the safe `name2rid` API and reused
  for every file, since the inputs share contig order. This gives index-seek performance
  without the memory overhead of a synced reader. **Region mode requires a `.csi` index.**

The merged header is the **first** input's header (template) plus the sample columns of
every subsequent input, appended in order. Records are rebuilt against that merged header;
`INFO`/`FORMAT` tag *types* are resolved once up front from the first input's header.

## 3. Command-line interface

```
paste-vcfs [--region chr:start-end[,chr2:...]] [--threads N] \
           [--info TAG1,TAG2,...] [--format TAG1,TAG2,...] \
           --output OUT.bcf  INPUT1 INPUT2 [INPUT3 ...]
```

| Flag / arg | Meaning | Default |
|---|---|---|
| `-r, --region` | One or more comma-separated regions to emit (index-seek per region). Coordinates are 1-based inclusive; commas inside coordinates are stripped, so `chr1:1,000-2,000` works. | whole file |
| `-t, --threads` | Writer compression threads (values `> 1` enable threading). | `1` |
| `--info` | Comma-separated `INFO` tags to copy from the **base** (first) file onto the output. | none |
| `--format` | Comma-separated `FORMAT` tags to paste from **every** file. `GT` is treated as integer-typed. | none |
| `-o, --output` | Output BCF path. | *required* |
| `INPUT…` | Two or more input BCF/`VCF.GZ` files (the first is the "base"). | *required* |

## 4. Output

A BCF with the union of all input samples (base samples first, then each side file's
samples in argument order). Per record it writes `CHROM/POS/ID/REF/ALT/QUAL` from the base
file, the requested `INFO` tags from the base file, and, for each requested `FORMAT` tag,
the concatenation of that tag's per-sample values across all inputs in argument order.

## 5. Mathematics and exact semantics

`paste-vcfs` performs **no arithmetic** — it is a structural/collation operation. The exact
semantics that matter for testing are:

### 5.1 Row-alignment assertions (per record, per side file)

The following mismatches each abort the run with a `FATAL` message:

* `rid`/`pos` mismatch → `Position mismatch!`
* `ID` mismatch → `ID mismatch!`
* REF/ALT (`alleles`) mismatch → `REF/ALT mismatch!`
* If a side file runs out of records before the base → `Side file ran out of records`.

### 5.2 `INFO` copy (base only)

For each `--info TAG`, the tag's value is read from the **base** record according to its
declared type (`Integer`/`Float`/`String`/`Flag`) and pushed onto the output only if
present. Flags are emitted only when set on the base.

### 5.3 `FORMAT` paste (all files, concatenated)

For each `--format TAG`, values are collected from the base file first, then from each side
file in order, concatenated into one buffer, and pushed onto the output. The tool enforces
**presence symmetry**: if the tag is present on the base it must be present on every side
file, and vice-versa; otherwise it aborts with a precise
`FORMAT tag '<TAG>' mismatch at chr:pos …` message naming the offending file. `GT` is
handled as an integer tag (encoded allele indices). Supported types are Integer, Float, and
String; any other type for a requested `FORMAT` tag aborts.

### 5.4 Region handling

Each `--region` entry is parsed into iterator bounds (half-open internally) and an
in-memory guard re-checks `rid`/`pos` so shard edges are exact even if the index
over-returns. Multiple regions are processed in the order given.

## 6. Notes for pipeline maintainers

* Inputs **must be pre-aligned**. `paste-vcfs` is intentionally *not* a merge tool: it will
  not reconcile differing site lists, it will abort. Produce the shards from a common site
  list.
* `--region` needs a **CSI** index (`bcftools index` / `bcftools index -c`), not TBI —
  the region reader hard-codes `HTS_FMT_CSI`.
* Only tags named on `--info` / `--format` are carried through; everything else on the
  input records is dropped. `QUAL` is taken from the base; `FILTER` is not propagated.
* The first input defines the output contigs, `INFO`/`FORMAT` header definitions, and the
  leading sample block.
