# FAN: Focal Amplification in Neochromosome


AmpliconClassifier 2.0 introduces a new mechanistic call, **FAN** (focal
amplification in neochromosome). This page explains what a FAN is, how it
differs from the amplification classes AC already reported, and what a FAN call
looks like in practice. It is intended as orientation for users encountering the
`FAN+` column for the first time. A full description of the phenotype and its
derivation is in preparation.

---

## 1. What is a FAN?

A **FAN** is a focal oncogene amplification that is **broad, structurally
complex, and multi-chromosomal**, and that resolves into a chromosomal or
**neochromosomal** structure rather than an extrachromosomal circle.

A *neochromosome* is a large, structurally rearranged chromosome assembled *de
novo* from amplified material, frequently drawing on several chromosomes of
origin and carrying amplified oncogenes (Garsed et al.,
[Cancer Cell 2014](https://doi.org/10.1016/j.ccr.2014.09.010)). FAN is the term
we introduce for the focal amplifications that sit inside such structures.

The defining contrast is with ecDNA. Both are focal oncogene amplifications, but
they are built differently:

| | **ecDNA**                                     | **FAN**                                                       |
|---|-----------------------------------------------|---------------------------------------------------------------|
| Physical form | Extrachromosomal circle                       | Chromosomally integrated / neochromosomal                     |
| Genomic territory | Compact, typically a few Mb                   | Broad, tens to hundreds of Mb                                 |
| Copy number | More likely to be uniform across the element  | **Heterogeneous** across the amplicon                         |
| Chromosomes of origin | Frequently one dominant locus                 | Two or more, each contributing substantial amplified material |
| SV content | Variable number of junctions relative to size | Dense, with many large or interchromosomal junctions          |

A useful intuition: **ecDNA tend to be sharper and better circumscribed; a FAN is wide and rugged.** The
single most common misreading of the FAN call is to assume it means "a very
large amplification." It does not. Copy number alone does not make a FAN; see
the ecDNA contrasts in §3.2, which include the *highest* copy number of any
example on this page and are firmly FAN-negative.

FAN is also **not mutually exclusive** with the other mechanistic calls. A FAN
event can incorporate ecDNA- and BFB-type mechanisms within the same
amplification, so `FAN+` can co-occur with `ecDNA+` and `BFB+` on one amplicon
(SNU16 amplicon1 is such a case in our cell-line panel; see §4). Treat `FAN+` as
an independent flag, not as a competing label.

---

## 2. How AmpliconClassifier calls FAN

Every amplicon is scored by a lightweight logistic-regression model that reads
the AmpliconArchitect breakpoint graph and computes five features:

| Feature | What it captures |
|---|---|
| `amplified_span_mb` | Genomic territory of the amplification: per-chromosome coordinate envelope of segments above the CN floor (default 4), summed over chromosomes. |
| `sv_qualifying_inter_or_large` | Count of SV junctions that are interchromosomal or intrachromosomal and >1 Mb (short <50 kb deletions/duplications are excluded). A measure of long-range SV burden. |
| `chrom_dominant_frac` | Fraction of amplicon sequence contributed by its single largest chromosome. **Low** values mean the material is spread across chromosomes. |
| `sv_crossing_frac` | Fraction of qualifying junctions that span at least one other junction's breakpoint, i.e. how interleaved, rather than nested and orderly, the rearrangements are. |
| `amplified_cn40_span_mb` | Envelope span of segments above CN 40. Acts mainly as a suppressor: compact, very-high-CN ecDNA concentrates its signal here. |

The model outputs a probability; `fan_call` applies a decision threshold of
**0.5**. All weights are compiled into `ampclasslib/fan_ml.py`. There is no
external model file and nothing to download or configure.

**How the model was built.** Labels came from manual review of AmpliconArchitect
reconstructions in existing pan-cancer focal-amplification calls
([Kim et al., *Nature Genetics*, 2024](https://doi.org/10.1038/s41588-024-01949-7))
and in cancer cell lines, with each amplicon annotated as FAN-like or not and
weighted by the annotator's confidence. Candidate models were fit and compared
under cross-validation grouped by sample, so that amplicons from the same tumor
never appeared in both training and evaluation, and were then scored on
held-out amplicons that played no part in model development. The feature set was
deliberately kept small in favor of interpretability. Orthogonal confirmation
comes from the cytogenetically characterized cell lines of Luebeck et al.,
*Nature Methods*, 2026 (in press), where a call can be checked against metaphase FISH; the
examples in §3 are drawn from that panel.

**Where the FAN results appear in AC outputs:**

- `[prefix]_amplicon_classification_profiles.tsv`: the `FAN+` column
  (`Positive` / `None detected`), alongside `ecDNA+` and `BFB+`.
- `[prefix]_fan_calls.tsv`: per-amplicon `fan_call`, `fan_probability`, and all
  five feature values. This is the file to consult when you want to understand
  *why* an amplicon scored the way it did.
- `[prefix]_result_table.tsv`: FAN-positive amplicons appear as a numbered
  feature `FAN_1`, with a `FAN probability` column.
- FAN intervals are written as `[sample]_amplicon[N]_FAN_1_intervals.bed`.

---

## 3. Examples

The examples below are cell lines from our validation panel, classified with
AmpliconClassifier 2.0.0 against hg38. Each row pairs the AmpliconArchitect
amplicon plot (copy number across the amplified segments, with SV junctions
drawn as arcs) with metaphase FISH from the same line.

In the FAN examples the probe
signal is **clustered and chromosome-associated** into large homogenously staining region (HSR) amplifications. In the ecDNA examples it is
**dispersed as many small, separate puncta** through the spread: double
minutes, physically independent of any chromosome.

### 3.1 Three FANs

![FAN examples: AmpliconArchitect plots and FISH for BT474, SJSA1, and SKBR3](images/fan_examples_aa_fish.png)

**BT474 amplicon6** (0.996) shows that high copy number is not required at all:
not one segment exceeds CN 40. chr17 and chr20 contribute almost equally, giving
the lowest `chrom_dominant_frac` in the panel, with *ERBB2* and *ZNF217*
amplified together across the two chromosomes.

**SJSA1 amplicon2** (probability 0.997): amplified material
across chr4, chr12, and chr13 with the densest SV burden in the panel, carrying
*MDM2*, *CDK4*, *DDIT3*, *GLI1*, and *HMGA2* together in one structure.

**SKBR3 amplicon1** (0.945) is the broadest, drawing on chr3, chr8, chr14, and
chr17 and co-amplifying *MYC* and *ERBB2*. Its peak copy number is *lower* than
SJSA1's; it is called on breadth and multi-chromosomal spread, not amplitude.



### 3.2 Two complex ecDNA, for contrast

![ecDNA contrast: AmpliconArchitect plots and MYC FISH for COLO320DM and SCLC21H](images/ecdna_contrast_aa_fish.png)

**COLO320DM amplicon4** (0.031) is canonical *MYC* ecDNA. Its copy number is the highest on this
page and three times SJSA1's peak, but the amplification
is compact, and nearly all of it sits above CN 40.

**SCLC21H amplicon1** (0.071) removes the easy explanation that FAN just means
"big": its 143 Mb envelope essentially matches SJSA1's 141 Mb, and it is highly
complex, yet it is firmly negative. `chrom_dominant_frac` is 0.84, meaning chr8
overwhelmingly dominates, against 0.49 for BT474's genuinely co-dominant pair.

Neither negative call means the amplicon is simple. Both of these are among the
most structurally complex amplicons in our collection of FISH-analyzed ecDNA.

### 3.3 The five amplicons side by side

| Amplicon | FAN prob. | span (Mb) | SVs | chrom_dominant | crossing | CN>40 span (Mb) |
|---|---|---|---|---|---|---|
| SJSA1 amplicon2 | 0.997 | 140.8 | 176 | 0.63 | 0.73 | 73.3 |
| SKBR3 amplicon1 | 0.945 | 248.5 | 68 | 0.78 | 0.64 | 12.5 |
| BT474 amplicon6 | 0.996 | 56.8 | 41 | 0.49 | 0.69 | 0.0 |
| COLO320DM amplicon4 | 0.031 | 3.7 | 12 | 0.61 | 0.62 | 3.4 |
| SCLC21H amplicon1 | 0.071 | 143.1 | 41 | 0.84 | 0.71 | 116.5 |

**A caution on reading span.** `amplified_span_mb` is a per-chromosome coordinate
envelope (§2), not a count of amplified bases. SCLC21H's 143 Mb envelope contains
only about **2 Mb** of sequence above the CN floor, scattered as widely separated
islands; SJSA1's near-identical envelope is filled by roughly **18 Mb** of real
amplified material, and BT474's much smaller envelope by **30 Mb**. A wide
envelope can mean a genuinely broad amplification or just a few distant
fragments; the other features are what tell them apart.

---

## 4. What changed from AmpliconClassifier 1.X

**All three FAN examples above were reported simply as `ecDNA+` by earlier
versions of AmpliconClassifier.** Under AC 1.5.X, SJSA1 amplicon2, SKBR3
amplicon1, and BT474 amplicon6 each carried `ecDNA+: Positive`; FAN did not
exist as a category, so a broad, cyclic-looking neochromosomal amplification had
nowhere else to go and was absorbed into the ecDNA call. In 2.0 the same three
amplicons are `FAN+: Positive` with `ecDNA+: None detected`.

The mechanism is a deliberately strict gate. When the FAN classifier fires, AC
no longer accepts the ordinary cyclic evidence as sufficient for an ecDNA call:
it retains only circular cycles larger than 250 kbp in which **every** segment's
graph copy number exceeds `fan_ecDNA_min_cn` (default **60**).

**ecDNA is very plausibly involved in forming a FAN**, since an extrachromosomal
intermediate that subsequently integrates is a natural route to a neochromosome,
and ecDNA may well persist in
rarer subclones alongside the neochromosomal structure even after the dominant
population has integrated it. What the reassignment says is narrower: the
structure best supported by *this* bulk reconstruction is neochromosomal rather
than a circle. A FAN call describes the dominant architecture the data support,
not the full evolutionary history of the amplification, and it should not be
read as excluding a past or minor extrachromosomal population. Single-cell or
cytogenetic data remain the way to detect a residual ecDNA subclone.

**If you are comparing cohorts across versions,** expect some amplicons to move
from `ecDNA+` to `FAN+`. This is the intended behavior change, not a
regression. Amplicons affected are broad and SV-dense, so the shift concentrates
in the complex tail of the cohort rather than in ordinary focal ecDNA.
`[prefix]_fan_calls.tsv` gives the per-amplicon detail needed to audit any
change you see.

---

## 5. Interpreting a FAN call

**Cytogenetic expectation.** In cell-line FISH validation, FAN calls
corresponded predominantly to non-native HSRs (amplified material integrated at
a chromosomal location other than the locus of origin), consistent with the
neochromosomal interpretation. No FAN call corresponded to a double minute
alone; where a FAN-classified amplicon coincided with a FISH-observed DM, that
amplicon was also called `ecDNA+`.

**A few practical notes:**

- `FAN+` is an independent mechanistic flag and can co-occur with `ecDNA+` and
  `BFB+`; do not treat the three as a single-choice label, and note that the
  abstract `amplicon_decomposition_class` is a separate axis again. Independent
  does not mean unrelated, though: as §4 describes, a FAN call does raise the
  evidence bar an ecDNA call must clear within that amplicon.
- `fan_probability` is a model score, not a calibrated posterior over
  biological ground truth. For borderline amplicons near 0.5, inspect the
  feature values in `_fan_calls.tsv` and the AA plot before committing to an
  interpretation.

---

## 6. Summary

FAN is a new category introduced with AmpliconClassifier 2.0 representing amplifications that appear embedded in neochromosomes. 
They are recognized by the combination of broad genomic territory, a dense burden
of large or interchromosomal rearrangements, and amplified material drawn from
more than one chromosome, a signature that separates them from the compact,
uniformly high-copy profile of ecDNA. For existing users the practical
consequence is that some amplicons previously reported as `ecDNA+` are now
reported as `FAN+`, with `[prefix]_fan_calls.tsv` giving the per-amplicon
features behind any individual call.
The full description of the phenotype, the model's derivation, and its validation are in
preparation; this page is intended only to give users a working sense of what
the `FAN+` flag means.

For the output-format details referenced above, see the
[README](../README.md#3-outputs).
