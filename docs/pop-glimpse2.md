# `pop-glimpse2`

Source: `resources/pop-glimpse2/src/bin/pop-glimpse2.rs`

## 1. High-level function

The imputation panel fundamentally consists of **multi-allelic "bubbles"**: each bubble site offers several
ALT haplotype *paths*, and each path is really a concatenation of one or more **atomic** bi-allelic
variants (SNPs/indels/SVs). See `docs/extract-bubble-PLs.md` for background on bubble representation.

For the purposes of using this panel with GLIMPSE2, we first split these multi-allelics to bi-allelics. This is
required by GLIMPSE2, which runs an HMM that is fundamentally bi-allelic. Such naive treatment of multi-allelics
when using GLIMPSE2 is standard. Unfortunately, this means that GLIMPSE2 can emit combinations of haplotypes that do not respect ploidy;
e.g., it could emit hom-alt for all paths in a given multi-allelic bubble. (Even though it would seem that the zero genetic
distance between all sites in a bubble should suppress such emission in the HMM, it seems in practice this is not the case;
perhaps the emission error term relaxes this constraint sufficiently.)

`pop-glimpse2` thus inverts the conversion to bubble representation, while simultaneously resolving consistent haplotypes
from all of the paths emitted by GLIMPSE2. Given GLIMPSE2's **phased
posterior probabilities** over the bubble paths, it computes, **for each haplotype
separately**, a probability distribution over the paths and then **redistributes that
probability onto the constituent atomic variants**, emitting one bi-allelic record per
atomic variant with a reconstructed `GT`, dosage `DS`, and genotype posteriors `GP`.

In one sentence: **it turns per-path phased GLIMPSE2 output back into per-atomic-variant
phased dosages and genotype probabilities.**

## 2. Design

The tool is a three-way **lockstep / windowed streaming** program:

* **stdin** — the multi-allelic imputed VCF (one line per path; paths that belong to the
  same bubble share a `POS`). Must carry `FORMAT/GT` and `FORMAT/GP` for every sample and
  `INFO` fields `RAF`, `AF`, `INFO`.
* **arg 1 `<biallelic ID VCF>`** — the atomic-variant dictionary. Each record's
  `INFO/ID=<atomic_id>` maps an atomic id to its true `(POS, ID, REF, ALT)`. Read into a
  windowed `id_buffer` (see `window_size`).
* **arg 2 `<sites VCF>`** — read **in exact lockstep** with stdin (line *i* of sites must
  have the same `CHROM/POS/REF/ALT` as line *i* of stdin). Its `INFO/ID=<a:b:c>` lists the
  atomic ids that make up each path. A mismatch aborts the run.

Because the three streams are coordinate-synchronised, the tool holds only one bubble's
worth of path lines plus a position-windowed slice of the atomic dictionary in memory.

`mimalloc` is installed as the global allocator for throughput.

Header lines from stdin are passed through unchanged **except** that `INFO=<ID=AK...`,
`FORMAT=<ID=GL...`, and `FORMAT=<ID=KC...` definitions are dropped (those tags are not
present on the output).

## 3. Command-line interface

```
cat <multiallelic VCF> | pop-glimpse2 <biallelic ID VCF> <sites VCF> \
    [max_alleles] [window_size]
```

| Position | Meaning | Default |
|---|---|---|
| stdin | Multi-allelic imputed VCF (paths), with `GT` + `GP` and `INFO` `RAF/AF/INFO`. | *required* |
| `<biallelic ID VCF>` (arg 1) | Atomic-variant dictionary: `INFO/ID` → `(POS, ID, REF, ALT)`. | *required* |
| `<sites VCF>` (arg 2) | Lockstep sites file whose `INFO/ID=<a:b:...>` lists each path's atomic ids. | *required* |
| `max_alleles` (arg 3) | Maximum number of paths kept **per haplotype per sample** (the top-scoring ones). | `10` |
| `window_size` (arg 4) | bp radius of the atomic-dictionary buffer around the current site. Must exceed the largest bubble span or an atomic id will be "not found". | `500000` |

## 4. Output

Plain VCF written to **stdout**: the retained header, then one bi-allelic record per atomic
variant (sorted by dictionary `POS`, then `REF`, then `ALT`) with
`FORMAT = GT:DS:GP` and `INFO = ID=<atomic_id>;RAF=…;AF=…;INFO=…` (the `RAF/AF/INFO` copied
from the first path that contains the atomic variant).

## 5. Mathematics

All arithmetic below is performed in **32-bit floating point** (`f32`); this matters for
bit-exact reproduction (see the testing notes).

### 5.1 Per-path phased haplotype probabilities

For sample *s* and path *a*, let GLIMPSE2's posteriors be `GP = (gp0, gp1, gp2)` =
`P(0/0), P(0|1 or 1|0), P(1/1)` and let the phased hard-call be the record's `GT`. Define
the probability that **hap 0** and **hap 1** carry this path's ALT as:

```
GT = "1|0":  (p0, p1) = (gp2 + gp1, gp2)
GT = "0|1":  (p0, p1) = (gp2,       gp2 + gp1)
otherwise :  (p0, p1) = (gp2 + gp1/2, gp2 + gp1/2)
```

Each is then clamped to `[1e-5, 1 − 1e-5]`. Intuitively `gp2` (hom-alt) contributes to both
haplotypes, while the heterozygous mass `gp1` is assigned to the phased haplotype (or split
evenly when phase is unknown).

### 5.2 Per-haplotype path selection and normalisation

For each sample and haplotype `h ∈ {0,1}` independently:

1. **Score** every path by `p0 + p1` and keep the top `m = min(#paths, max_alleles)`.
2. Convert each kept path's probability to **odds** `w_h(a) = p_h(a) / (1 − p_h(a))`.
3. Normalise with an **implicit reference pseudo-path of weight 1**:

```
Z_h = 1 + Σ_a w_h(a)
p̂_h(a) = w_h(a) / Z_h
```

This produces a proper multinomial distribution over `{REF, kept paths}` on each
haplotype: the "1" reserves probability for "no ALT path taken", and the odds transform
makes the selection scale-free.

### 5.3 Projection onto atomic variants

Each atomic variant `v` receives, per haplotype, the summed normalised probability of every
kept path that contains it:

```
P_h(v) = Σ_{a : v ∈ path a} p̂_h(a)
```

accumulated across the (up to `max_alleles`) kept paths, per sample. Call the two results
`p0 = P_0(v)` and `p1 = P_1(v)`, each finally clamped to `[0, 1]`.

### 5.4 Emitted `GT`, `DS`, `GP`

Assuming the two haplotypes are independent Bernoulli(`p0`), Bernoulli(`p1`):

```
GT   : hap0 = 1 if p0 > 0.5 else 0 ; hap1 = 1 if p1 > 0.5 else 0   (phased, "h0|h1")
DS   : p0 + p1                                                      (expected ALT count)
GP0  : (1 − p0)(1 − p1)          # P(0/0)
GP1  : p0(1 − p1) + (1 − p0)p1   # P(0/1)
GP2  : p0·p1                     # P(1/1)
```

The three `GP` values are rendered as integers per mille:
`v_i = round(GPi · 1000)` (round half away from zero). Any rounding residual
`1000 − Σ v_i` is added to the **largest** bin so the trio sums to exactly 1000, then each
is printed as `v_i / 1000`.

### 5.5 Number formatting

`format_float(x)` prints `x` with three decimals then strips trailing zeros and a trailing
`.`, mapping the empty result to `"0"` (so `0.500 → "0.5"`, `1.000 → "1"`, `0.000 → "0"`).
`DS` is `format_float(p0 + p1)`; each `GP` component is `format_float(v_i / 1000)`.

## 6. Notes for pipeline maintainers

* The **lockstep invariant** between stdin and the sites VCF is strict: any
  `CHROM/POS/REF/ALT` divergence, or the sites file running short, aborts with a `FATAL`
  message. Generate both files from the same source pass.
* Every atomic id listed in a path's `INFO/ID` **must** be resolvable in the id buffer;
  if `window_size` is smaller than the distance between a bubble and its atomic
  dictionary entries the run aborts with "Variant ID … not found in the ID buffer".
* Because the numeric output is `f32`, treat it as a **characterization ("golden")**
  target for exact-match testing: freeze a known-good output and diff against it. A
  re-implementation in another language's default `f64`/`round-half-to-even` will differ in
  the last digit on ties. The test harness in this repo does exactly this (see
  `resources/tests/pop-glimpse2/`).
