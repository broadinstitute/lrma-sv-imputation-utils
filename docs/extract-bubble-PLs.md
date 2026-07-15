# `extract_bubble_PLs`

Source: `resources/extract-bubble-PLs/src/main.rs`

## 1. High-level function

`extract_bubble_PLs` builds the **per-sample likelihood evidence** that GLIMPSE2 needs
in order to impute a cohort against a reference panel.

The reference panel is a **bi-allelic BCF** in which every record is one "bubble" — a
single REF/ALT site that GLIMPSE2 will impute. For each panel bubble the tool walks an
input callset (either a per-sample **gVCF** or a multi-sample **joint VCF**) and copies
the genotype likelihoods (`PL`) for that exact allele onto the panel coordinate, emitting
a new bi-allelic BCF whose sites are *identical to the panel* but whose `FORMAT` carries
`GT` and a rescaled `PL` for every requested sample.

In one sentence: **it re-expresses the input callset's likelihoods on the fixed panel
coordinate grid, one bi-allelic row per panel site.**

The output BCF is the direct input to GLIMPSE2's likelihood phase.

## 2. Design

The tool is a **streaming, index-seeking sweep** with a bounded look-ahead buffer, so it
never loads a whole chromosome into memory:

1. Open the panel with an `IndexedReader` (it can be sharded by `--region`).
2. Open the input (gVCF/joint VCF) with an `IndexedReader`.
3. Build the output header from the **panel** header (so contigs/IDs come from the panel),
   then append `GT` and `PL` `FORMAT` lines and the requested sample columns.
4. Iterate panel records in order. A `VecDeque<InputRec>` holds a sliding window of input
   records around the current panel position (`--window` bp, plus a fixed 50 bp slack and,
   for gVCF mode, a 100 kb backward padding to catch upstream reference blocks / shifted
   indels). The buffer is *advanced* (read more input) and *pruned* (drop records that
   fall behind the panel cursor) on every panel site.
5. For each panel bubble, search the buffer for a matching ALT and emit `GT`/`PL`; if
   nothing matches, apply the fallback rule for the mode (below).

Contig identity between the two files is resolved through a `chrom_rank` map built from the
**input** header, so the two coordinate systems can be compared even when their contig
orderings differ.

### Custom index adapter (`##idx##`)

Both the panel and input path arguments accept a Cromwell/WDL-style
`path/to/file.bcf##idx##path/to/file.bcf.csi` suffix. When present, the tool symlinks the
provided index next to the data file (if the expected `.csi`/`.tbi` name does not already
exist) so that HTSlib's index auto-detection succeeds. This lets the tool run against
localized files whose indexes were delocalized to a different name.

## 3. Command-line interface

```
extract_bubble_PLs <gvcf|joint> <panel.bcf> <input.vcf.gz> <output.bcf> \
    [--region chr:start-end] [--samples list.txt] [--window 15000] \
    [--cap-pl 30] [--scale-pl 5.0] [--threads 4]
```

| Position / flag | Meaning | Default |
|---|---|---|
| `mode` (positional 1) | `gvcf` or `joint` — **must** be exactly one of these. Controls padding and the no-match fallback (see §5). | *required* |
| `panel.bcf` (positional 2) | Bi-allelic, **indexed** panel BCF. Supplies output contigs, coordinates, IDs, and REF/ALT. May carry the `##idx##` suffix. | *required* |
| `input.vcf.gz` (positional 3) | The **indexed** input callset (gVCF or joint VCF) to harvest likelihoods from. May carry the `##idx##` suffix. | *required* |
| `output.bcf` (positional 4) | Path of the BCF to write. | *required* |
| `--region chr:start-end` | Restrict processing to one contig / interval (1-based, inclusive). `chr` alone, `chr:pos`, and `chr:start-end` are all accepted. Used for sharding. | whole file |
| `--samples list.txt` | Newline-delimited list of sample names to keep (subset + reorder of the input samples). Unknown names are warned and skipped. | all input samples |
| `--window` | Half-width, in bp, of the input look-ahead/look-behind buffer around the current panel site. | `15000` |
| `--cap-pl` | Upper clamp applied to every emitted `PL` value after scaling. | `30` |
| `--scale-pl` | Divisor applied to every input `PL` before capping. | `5.0` |
| `--threads` | HTSlib (de)compression threads for the readers and the writer. | `4` |

## 4. Output

A bi-allelic BCF with exactly the panel's contigs and one record per **processed** panel
site (see §5 for which sites are dropped). Each record carries:

* `ID` — copied from the panel record's `INFO/ID` string if present, else `.`.
* `FORMAT/GT` — a diploid hard-call, one of `0/0`, `0/1`, `./.` (always unphased).
* `FORMAT/PL` — three values `[PL(0/0), PL(0/1), PL(1/1)]`, each rescaled/capped.

## 5. Mathematics and exact semantics

### 5.1 Minimal representation

Before comparing an input ALT to a panel ALT, both are reduced to a canonical minimal
representation `get_minimal_representation(pos, ref, alt)`:

1. If the ALT begins with `<` (symbolic, e.g. `<NON_REF>`), return unchanged.
2. **Right-trim:** while REF and ALT share a common last base, drop it from both.
3. **Left-trim:** while REF and ALT share a common first base, drop it from both and
   advance `pos` by 1.

Two variants "match" iff their minimal `(pos, ref, alt)` triples are byte-for-byte equal.

### 5.2 PL genotype indexing

A multi-allelic input record with `n` alleles stores `PL` (or `LPL`) as a flat vector of
length `n(n+1)/2`, in the VCF "GL ordering" where the genotype with alleles `a ≤ b` lives
at index

```
idx(a, b) = b·(b + 1) / 2 + a
```

To lift the panel's bi-allelic bubble whose ALT equals input allele index `k`
(`k = k_idx + 1`, the 1-based ALT position among the input alleles), the tool reads three
entries:

```
idx0 = idx(0, 0) = 0            # hom-ref
idx1 = idx(0, k) = k·(k + 1)/2  # het
idx2 = idx(k, k) = idx1 + k     # hom-alt
```

These become the output `PL[0..3]`. This is why `pl_stride = n(n+1)/2` is recomputed per
record and the read is guarded by `pl_stride > idx2`.

The genotype is remapped the same way: an input allele of `0` maps to output allele `0`,
an input allele equal to `k` maps to output allele `1`, and anything else maps to missing.

### 5.3 PL rescaling

Every emitted `PL` value passes through:

```
process_pl(v) = v                      if v == i32::MIN          (missing sentinel, preserved)
                min(cap_pl, trunc(v / scale_pl))   otherwise
```

`trunc(v / scale_pl)` is integer truncation toward zero of the floating-point quotient
(`((v as f64) / scale_pl) as i32`). With the defaults `scale_pl = 5.0`, `cap_pl = 30`, an
input `PL` of `99` becomes `min(30, trunc(19.8)) = 19`, and `255` becomes `min(30, 51) = 30`.

Rescaling compresses the dynamic range of the likelihoods (GLIMPSE2 does not need the full
Phred range and smaller values speed up its arithmetic) and the cap bounds the maximum
penalty any single site can contribute.

### 5.4 Match and fallback logic per panel site

For each panel bubble the output `GT`/`PL` are first reset to *missing / all-zero*, then:

1. **Variant match.** Scan the buffered input records on the same contig. For each
   non-reference-block record, compare each ALT's minimal representation to the panel
   bubble's. On the first equal match, copy `GT` and the three `PL` entries (via §5.2/§5.3)
   for every output sample, and stop.

2. **No match → fallback depends on mode:**
   * **`joint`** — the row is **dropped entirely** (`continue`). A bubble absent from the
     joint callset contributes nothing, which saves disk and GLIMPSE2 compute.
   * **`gvcf`** — search the buffer for a reference block (`<NON_REF>` with an `END`) that
     spans the site. If found, emit a confident hom-ref call: `GT = 0/0` and
     `PL = [0, process_pl(GQ), process_pl(2·GQ)]`, i.e. the reference-block `GQ` is used as
     the het penalty and twice `GQ` as the hom-alt penalty. If no covering block is found,
     the site is emitted as missing (`GT = ./.`, `PL = 0,0,0`).

So in **`gvcf`** mode every panel site produces a row (call, hom-ref, or missing); in
**`joint`** mode only sites present in the joint callset produce a row.

### 5.5 Region filtering

With `--region`, the panel reader is `fetch`ed to the interval and an additional
`pos < start || pos > end` guard drops any records the index over-returns, so shard
boundaries are exact.

## 6. Notes for pipeline maintainers

* The panel and input **must be indexed** (`.csi` or `.tbi`); the tool opens both with an
  `IndexedReader`.
* `PL` is preferred from the `LPL` (local-allele PL) `FORMAT` tag if present, otherwise the
  standard `PL` tag.
* Output allele set per row is exactly `[panel REF, panel ALT]`; multi-allelic panels are
  not supported (records with `< 2` alleles are skipped).
* The tool writes progress and shard diagnostics to `stderr` (prefixed `[DEBUG]`,
  `[INFO]`, `[Progress]`); only the final `[INFO] Preprocessed imputation BCF generated
  successfully.` line goes to `stdout`.
