# AmpliconClassifier Configuration Parameters

This document describes all configuration parameters available in AmpliconClassifier's `default_config.json` file. These parameters control various aspects of amplicon classification, including size thresholds, copy number cutoffs, and structural variant detection criteria.

---

## Table of Contents
- [General Amplicon Parameters](#general-amplicon-parameters)
- [Copy Number Parameters](#copy-number-parameters)
- [Cycle Classification Parameters](#cycle-classification-parameters)
- [Segmental Duplication Parameters](#segmental-duplication-parameters)
- [Decomposition Parameters](#decomposition-parameters)
- [BFB (Breakage-Fusion-Bridge) Parameters](#bfb-breakage-fusion-bridge-parameters)

---

## General Amplicon Parameters

### `tot_min_del`
- **Default:** 5000 (bp)
- **Description:** Minimum size threshold for a deletion to be considered non-trivial or significant. Used when determining if a cycle is rearranged by checking distances between consecutive segments.
- **Usage:** In `isRearranged()`, if the distance between two consecutive segments exceeds this threshold, the cycle is classified as rearranged.

### `minCycleSize`
- **Default:** 5000 (bp)
- **Description:** Minimum size of an AmpliconArchitect (AA) cycle to be considered as a valid amplicon. Cycles smaller than this are typically filtered out as too small to represent true amplifications.
- **Usage:** 
  - Used throughout the code to filter small cycles
  - In `cycleIsNoAmpInvalid()` to reject cycles below this threshold
  - Can be overridden with `--min_size` command-line argument

---

## Copy Number Parameters

### `min_amp_cn`
- **Default:** 4.5
- **Description:** Minimum copy number threshold for a segment to be considered part of a focal amplification. This is the baseline CN requirement for calling an amplicon.
- **Usage:** Used in multiple functions including `cycleIsNoAmpInvalid()` and `classifyBFB()` to determine if copy numbers are sufficiently elevated to indicate amplification.

### `sig_amp`
- **Default:** 7
- **Description:** Minimum copy number for a segment to be considered significantly amplified. This is a higher threshold than `min_amp_cn` and indicates stronger amplification.
- **Usage:** 
  - Used in `cycleIsNoAmpInvalid()` to determine if cycles show sufficient CN elevation
  - Used in `segdup_cycle()` to exclude high-CN regions from being classified as segmental duplications

### `high_amp`
- **Default:** 12
- **Description:** Minimum copy number threshold for high-level amplification. This represents very strong amplification events.
- **Usage:** Used in `segdup_cycle()` to exclude extremely high CN regions from being called as segmental duplications.

### `ampLenOverMinCN`
- **Default:** 5000 (bp)
- **Description:** Minimum amount (length) of amplicon sequence that must be above `min_amp_cn` to trigger a focal amplification call. Ensures amplified regions are of sufficient size to be biologically meaningful.
- **Usage:** Used in `classifyAmpliconProfile()` and `classifyBFB()` to validate that amplicons have sufficient high-CN content.

---

## Cycle Classification Parameters

### `compCycContCut`
- **Default:** 50000 (bp)
- **Description:** Minimum total size of complex cycles for classifying an amplicon as containing complex cyclic structures. This is specifically for complex (rearranged) cycles.
- **Usage:** Used in `classifyAmpliconProfile()` to determine if complex cyclic content is substantial enough to classify the amplicon as "Cyclic".

### `anyCycContcut`
- **Default:** 10000 (bp)
- **Description:** Minimum total size of any cyclic path/cycle (simple or complex) for an amplicon to be considered as containing meaningful cyclic content.
- **Usage:** Used in `classifyAmpliconProfile()` to ensure cyclic amplicons have sufficient cyclic content.

### `cycCut`
- **Default:** 0.12
- **Description:** Minimum proportion (fraction) of cycle weight in the decomposition for classifying an amplicon as ecDNA. This is the threshold for the combined weight of trivial and complex cycles relative to total amplicon weight.
- **Usage:** In `classifyAmpliconProfile()`, if `cycSig` (sum of "Trivial cycle" and "Complex-cyclic" weights) exceeds this threshold (along with size requirements), the amplicon is classified as "Cyclic".

### `compCut`
- **Default:** 0.3
- **Description:** Minimum proportion of cycle weight attributed to "complex" (complex non-cyclic and cyclic) structures for classification purposes.
- **Usage:** Used in `classifyAmpliconProfile()` as an alternative pathway for classification when cyclic signatures don't dominate.

---

## Segmental Duplication Parameters

### `max_segdup_size`
- **Default:** 1000000 (bp / 1 Mb)
- **Description:** Maximum allowed size for a cycle to potentially be classified as a segmental duplication. Cycles larger than this cannot be segmental duplications regardless of other properties.
- **Usage:** In `segdup_cycle()`, cycles exceeding this size are immediately excluded from being called as segmental duplications.

### `segdup_max_extra_fraction`
- **Default:** 0.25
- **Description:** Maximum additional copy number ratio beyond the baseline 2× (diploid) for a region to be classified as a segmental duplication. For example, with default value 0.25, the cycle CN / border CN ratio must be ≤ 2.25 to be called a segmental duplication.
- **Usage:** In `segdup_cycle()`, the ratio of cycle CN to flanking CN is calculated. If `(ratio - 2) <= segdup_max_extra_fraction`, and other criteria are met, the cycle may be called a linear segmental duplication rather than an ecDNA amplification.

---

## Decomposition Parameters

### `decomposition_strictness`
- **Default:** 0.1
- **Description:** Scaling factor (between 0 and 1) controlling how strictly to filter low copy number decomposed paths/cycles. Higher values filter more aggressively, removing more low-weight decompositions. When multiplied against the maximum copy number found in the amplicon, it defines the minimum assigned CN a trivial valid cycle can have. Used to control the CN threshold for singleton paths/cycles when determining validity.
- **Usage:** 
  - In `cycleIsNoAmpInvalid()`, multiplied by `maxCN` to create a threshold: `scale = min(min_flow, maxCN * decomposition_strictness)`
  - Can be overridden with `--decomposition_strictness` command-line argument
  - Default can be changed via config file

### `min_flow`
- **Default:** 1.0
- **Description:** Minimum flow (copy number weight) to consider a path in the AA decomposition as valid for classification as a focal amplification. Paths with flow below this threshold are filtered out.
- **Usage:** 
  - Used in `cycleIsNoAmpInvalid()` as part of the scaling calculation
  - Used in `clusterECCycles()` to filter cycles with insufficient CN
  - Can be overridden with `--min_flow` command-line argument

---

## BFB (Breakage-Fusion-Bridge) Parameters

BFB is a specific type of genomic amplification mechanism characterized by inverted duplications and foldback structures.

### `min_fb_read_prop`
- **Default:** 0.25
- **Description:** Minimum proportion of structural variant reads supporting foldback inversions required to call a BFB event. This ensures sufficient read-level evidence for the foldback structures characteristic of BFB.
- **Usage:** In `classifyBFB()`, if `fb_read_prop < min_fb_read_prop`, BFB classification is rejected.

### `fb_break_weight_prop`
- **Default:** 0.3
- **Description:** Minimum proportion of utilized path/cycle structural variants that must support BFB foldback patterns, where each SV is weighted by the path/cycle's assigned CN. This measures how dominant the BFB signature is among all structural variants in the amplicon, weighted by how dominant the path/cycle is.
- **Usage:** In `classifyBFB()`, if `fb_bwp < fb_break_weight_prop`, BFB classification is rejected.

### `fb_dist_cut`
- **Default:** 25000 (bp)
- **Description:** Maximum distance between breakpoint ends for an inversion to be classified as a foldback. Foldbacks are short-range inversions, so inversions spanning longer distances are not considered foldbacks.
- **Usage:** 
  - In `compute_f_from_AA_graph()`, inversions are only counted as foldbacks if they are on the same chromosome and within this distance
  - In `cycles_file_bfb_props()`, used to distinguish foldback breaks from distal breaks

### `max_nonbfb_break_weight`
- **Default:** 0.5
- **Description:** Maximum proportion of utilized path/cycle structural variants that can support non-foldback SV connections while still calling BFB, where each SV is weighted by the path/cycle's assigned CN. If too many SVs are non-foldback, the amplicon likely isn't a pure BFB.
- **Usage:** In `classifyBFB()` (via `cycles_file_bfb_props()`), if `nonbfb_sig > max_nonbfb_break_weight`, additional checks on cycle ratio are required.

### `min_bfb_cycle_weight_ratio`
- **Default:** 0.6
- **Description:** Minimum cycle CN × length-weighted proportion of BFB-like paths/cycles required for BFB classification when non-BFB content is also present. This ensures BFB structures dominate the amplicon.
- **Usage:** In `classifyBFB()`, when non-BFB signal is high, `bfb_cyc_ratio` must exceed this threshold to maintain BFB classification.

---

## Customizing Parameters

### Via Config File
You can create a custom configuration file and pass it using:
```bash
amplicon_classifier.py --config my_custom_config.json ...
```

### Via Command Line
Some parameters can be overridden directly via command-line arguments:
- `--min_flow`: Override `min_flow`
- `--min_size`: Override `minCycleSize`
- `--decomposition_strictness`: Override `decomposition_strictness`

### Example Custom Config
```json
{
    "minCycleSize": 10000,
    "min_amp_cn": 5.0,
    "sig_amp": 8,
    "cycCut": 0.15,
    "decomposition_strictness": 0.15
}
```

---

## Parameter Interactions

### Amplicon Size Classification
Multiple parameters work together to define amplicon size thresholds:
- `minCycleSize`: Base minimum size
- `anyCycContcut`: Minimum cyclic content
- `compCycContCut`: Minimum complex cyclic content
- `ampLenOverMinCN`: Minimum high-CN content

### Copy Number Hierarchy
Three CN thresholds create a hierarchy of amplification strength:
1. `min_amp_cn` (4.5): Basic amplification
2. `sig_amp` (7): Significant amplification
3. `high_amp` (12): High-level amplification

### BFB Detection
BFB classification requires multiple criteria to be met simultaneously:
1. Sufficient foldback read support (`min_fb_read_prop`)
2. Dominant foldback breakpoint weight (`fb_break_weight_prop`)
3. Limited non-foldback content (`max_nonbfb_break_weight`)
4. Strong BFB cycle representation (`min_bfb_cycle_weight_ratio`)
5. Adequate amplification (`min_amp_cn`, `ampLenOverMinCN`)

---

## Notes

- All size parameters are in base pairs (bp)
- All proportion parameters are fractions between 0 and 1
- Copy number parameters are in absolute copy number units
- Parameters marked with command-line override options provide runtime flexibility
- Most parameters have been empirically tuned based on analysis of real cancer genomics data

