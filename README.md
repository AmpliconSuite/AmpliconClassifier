# AmpliconClassifier
![GitHub release (latest by date)](https://img.shields.io/github/v/release/AmpliconSuite/AmpliconClassifier?display_name=release)
![GitHub](https://img.shields.io/github/license/AmpliconSuite/AmpliconClassifier)

### Classify [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) outputs to predict types of focal amplifications present.

This tool classifies the outputs of [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconSuite-pipeline).

If using AmpliconClassifier, please cite the following publication which describes the AmpliconClassifier methodology in the [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05937-5/MediaObjects/41586_2023_5937_MOESM1_ESM.pdf) section:

Luebeck J, et al. [Extrachromosomal DNA in the cancerous transformation of Barrett’s oesophagus](https://www.nature.com/articles/s41586-023-05937-5). Nature. 2023. PMID: 37046089

See the section at the end of the README for information about the legacy version.<br />
<br />

### Amplification types detected by AmpliconClassifier

AmpliconClassifier characterizes **focal oncogene amplifications** — the local, high-copy-number amplification of a genomic region — as reconstructed by AmpliconArchitect or CoRAL. It does **not** analyze karyotypic abnormalities such as aneuploidy, balanced translocations, or broad chromosome/arm-level gains and losses; it focuses on the structure of the focal amplification event itself.

The main amplification mechanisms and structures it distinguishes are:

- **ecDNA (extrachromosomal DNA)** — Circular, extrachromosomal amplicons that replicate autonomously and segregate unequally to daughter cells, driving rapid copy-number escalation and intratumoral heterogeneity. ecDNA is typically focal and high copy-number, and is reconstructed as one or more genome cycles. Reported with the `ecDNA+` flag.
- **BFB (breakage–fusion–bridge)** — Amplification produced by iterative cycles of chromosome-end breakage, sister-chromatid fusion, and anaphase-bridge rupture. BFB leaves a characteristic signature of fold-back inversions and a stair-step copy-number gradient, and the amplified material remains chromosomally integrated (frequently as a homogeneously staining region, HSR). Reported with the `BFB+` flag.
- **FAN (focal amplification in neochromosome)** — A broad, structurally complex, multi-chromosomal amplification that is integrated into a chromosome or assembled into a neochromosome rather than existing as an extrachromosomal circle. FAN is a term we introduce to describe these neochromosome-like focal amplifications. FAN events can incorporate ecDNA- and BFB-type mechanisms within the same amplification; `FAN+` is an independent flag and can co-occur with `ecDNA+` and `BFB+` (see the `_fan_calls.tsv` output below).
- **CNC (complex non-cyclic)** — An abstract label for a focal amplification embedded in a complex genome rearrangement that lacks the genome cycles characteristic of ecDNA. It describes the observed structure rather than a specific mechanism — the underlying process (e.g. chromothripsis) is often not fully clear.
- **Linear** — A focal amplification with few or no significant rearrangements, including low-copy amplifications arising from tandem duplication. The underlying mechanism is frequently unclear.

When AmpliconArchitect or CoRAL is run against the `GRCh38_viral` reference, an additional **Virus** category applies: amplicons corresponding to a viral genome are labeled `Virus`, and amplicons that fuse viral and human sequence are treated as human–viral **hybrids**. A hybrid oncoviral amplicon is reported through its human amplification mechanism (e.g. `ecDNA`) rather than as `Virus`, so its viral content is surfaced by the `contains_viral` column in the classification profiles and results table (`True` whenever any viral sequence is present, including hybrids).

**Mechanistic calls vs. the abstract decomposition class.** AmpliconClassifier reports these results in two complementary ways. The `ecDNA+`, `BFB+`, and `FAN+` columns — together with `Virus` — are *mechanistic* calls: each fires when evidence for that specific mechanism or origin is present, and more than one can be positive on the same amplicon. Separately, every amplicon also receives a single, abstract `amplicon_decomposition_class` that summarizes its overall structure. `Cyclic` is the structural signature underlying ecDNA and BFB calls, while `Linear` and `Complex-non-cyclic` are *residual* structural labels — they are what remains to describe an amplicon's structure (simple versus rearranged) when no mechanistic call fires, and they do not themselves imply a mechanism. The full set of `amplicon_decomposition_class` values is listed under [Outputs](#3-outputs) below.

### 1. Installation

AmpliconClassifier is included with [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline), but for standalone re-classification you can install it directly using the steps below.

#### Step 1: Create and activate a conda environment
```bash
conda create -n ampliconclassifier "python>=3.8"
conda activate ampliconclassifier
```

#### Step 2: Install dependencies
```bash
# Required
conda install -c conda-forge -c bioconda intervaltree scipy pandas numpy

# Optional: needed only for check_SV_support.py and BAM-based SV validation
conda install -c bioconda pysam

# Optional: needed only for coordinate lifting utility scripts
conda install -c conda-forge pyliftover
```

[BFBArchitect](https://github.com/AmpliconSuite/BFBArchitect), used for additional BFB detection, is a core dependency. Its published package is available on PyPI (not conda-forge/bioconda), so the standard installation in Step 3 resolves it with pip. Developers may instead install a local checkout in editable mode.

If you prefer pip for everything, clone the repository first and then install from the source tree in Step 3. The pip dependency list is also available in `requirements.txt`:

```bash
pip install -r requirements.txt
# optional: pip install pysam pyliftover
```

#### Step 3: Clone and install AmpliconClassifier
```bash
git clone https://github.com/AmpliconSuite/AmpliconClassifier.git
cd AmpliconClassifier
python -m pip install -e .
```

This normally installs BFBArchitect from PyPI, since it is a core dependency; a developer environment may substitute an editable local BFBArchitect checkout. BFBArchitect-positive regions are integrated into `BFB+` calls and reported BFB feature intervals; use `--no_bfbarchitect` to disable this integration at runtime. AC records the detected BFBArchitect version in its log. Passing reconstructions always write their graph, cycles, and visualization files under `bfbarchitect_outputs/`. The default classification profile reports the best BFBArchitect score and BFB call source; `--verbose_classification` adds raw classifier scores and extended BFBArchitect diagnostics.

Supported AmpliconSuite installation and container paths wire BFBArchitect's **Gurobi**, **Mosek**, and free, open-source **CBC** backends. Users choose which commercial capability to enable by which license, if any, they provide: Gurobi is recommended for large cohorts (roughly 50 or more samples) because it was fastest in our tests; Mosek is a supported alternative; and CBC requires no license. BFBArchitect automatically tries Gurobi, then Mosek, then CBC, warning and retrying if a selected commercial solver later fails. Gurobi autodetection verifies that its environment starts, but a restricted license may still reject a larger model at solve time and trigger fallback. Standalone AC developers must install the Mosek Python package separately if they want that optional backend. In our testing the BFB calls agree across solvers in the large majority of cases, though in rare borderline instances the choice of solver can shift a BFB score across the `BFB+`/`BFB-` threshold. See the [AmpliconSuite-pipeline README](https://github.com/AmpliconSuite/AmpliconSuite-pipeline#optimizer-licenses--do-i-need-one-short-answer-no) for details.

For older workflows that call scripts directly from the source tree, you may still set `$AC_SRC`:

```bash
echo export AC_SRC=$PWD >> ~/.bashrc
source ~/.bashrc
```

#### Step 4: Set up the AA data repo
AmpliconClassifier reads a small set of reference annotation files (mappability/low-complexity regions, gene models, oncogene list, centromeres) from a per-genome-build directory pointed to by `$AA_DATA_REPO`. First, set the environment variable and create the directory if needed:

```bash
echo 'export AA_DATA_REPO=~/data_repo/' >> ~/.bashrc
source ~/.bashrc
mkdir -p $AA_DATA_REPO
```

**Easiest: if you have [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) installed**, fetch a reference build with:
```bash
AmpliconSuite-pipeline.py --download_repo GRCh38   # or hg19, GRCh37, mm10, GRCh38_viral
```

**Without AmpliconSuite-pipeline**, the same data is a plain, unauthenticated download — no need to install AmpliconSuite-pipeline just for this:
```bash
ref=GRCh38   # or hg19, GRCh37, mm10, GRCh38_viral
wget https://refs.ampliconrepository.org/data/module_support_files/AmpliconArchitect/${ref}.tar.gz -P $AA_DATA_REPO
tar -xzf $AA_DATA_REPO/${ref}.tar.gz -C $AA_DATA_REPO
rm $AA_DATA_REPO/${ref}.tar.gz
```

Use the plain reference name (e.g. `GRCh38`), not the `_indexed` variant — the `_indexed` builds add a BWA index used only for AmpliconArchitect's read alignment step, which AmpliconClassifier never needs.

Mac users will also need:
```bash
brew install coreutils
```


### 2. Usage

`amplicon_classifier.py` takes a collection of (or single) AA graph files and corresponding AA cycles file as inputs.

**Most common - classifying multiple amplicons:**
You can provide the directory containing multiple AA amplicons or multiple uniquely named samples
>`python amplicon_classifier.py --ref GRCh38 --AA_results /path/to/AA/output/directories/ > classifier_stdout.log`

AC will crawl the given location and find all relevant AA files and perform classification on them.
Amplicons can be classified in parallel with `--jobs N`; output order follows the input file order.
When BFBArchitect is enabled, the approximate BFBArchitect solver thread budget is
`--jobs * --bfb_threads`. For multi-amplicon runs, `--bfb_threads 2` is usually
the most efficient per-amplicon setting; prefer increasing `--jobs` rather than
assigning many solver threads to each individual amplicon.

**To classify a single amplicon**:

>`python amplicon_classifier.py --ref GRCh38 --cycles sample_amplicon1_cycles.txt --graph sample_amplicon1_graph.txt > classifier_stdout.log`


**Less common - separate usage of `make_input.sh`:**

Alternatively, you can use the `make_input.sh` script to gather the necessary input files outside of AC:

`make_input.sh` takes a path and an output prefix. e.g:

>`make_input.sh /path/to/AA/output/directories/ example_collection` 

This would create a file called `example_collection.input` which can be given as the `--input` argument for AC.
The script pairs each `*_cycles.txt` file with the exact sibling `*_graph.txt` file sharing the same prefix,
reports missing graph files as errors, and reports orphan graph files as warnings.


#### Combining classification results from GRCh37 and hg19:

If combining data from both GRCh37 and hg19 in the same classification run, you can set the flag `--add_chr_tag` to add the "chr" prefix to each chromosome name and effectively unify everything as hg19-based.

### 3. Outputs

#### ****`[prefix]_amplicon_classification_profiles.tsv`**** 

Contains an abstract classification of the amplicon, and also indicates in separate columns "BFB+", "ecDNA+", and "FAN+" (Focal amplification in neochromosome) status.
Note that amplicons receiving a "Cyclic" classification may be ecDNA+, BFB+ or both.

| Column name                    | Contents                                                                                                                                                                                |
|--------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample_name`                  | Sample name prefix                                                                                                                                                                      |
| `amplicon_number`              | AA amplicon index, e.g. `amplicon2`                                                                                                                                        |
| `amplicon_decomposition_class` | Abstract description of the AA amplicon type.                                                                                   |
| `ecDNA+`                       | Prediction about whether the AA amplicon contains ecDNA. Note, an AA amplicon may contain regions surrounding the ecDNA, or multiple linked ecDNA. Either `Positive` or `None detected` |
| `BFB+`                         | Prediction about whether the AA amplicon is the result of a BFB. Either `Positive` or `None detected`                                                                                   |
| `FAN+`                          | Prediction about whether the AA amplicon is a FAN (Focal amplification in neochromosome): a broad, SV-dense, multi-chromosomal amplification that resolves into chromosomal/neochromosomal (non-circular) structure, as opposed to circular ecDNA. Reported as an independent flag and **not mutually exclusive** with `ecDNA+` or `BFB+`. Either `Positive` or `None detected`. See the `_fan_calls.tsv` description below. |
| `ecDNA_amplicons`              | Predicted number of distinct (non-overlapping) ecDNA which are represented in a single AA amplicon. This estimate is experimental.                                               |
| `contains_viral`               | Whether the amplicon contains any oncoviral sequence (a segment on a viral contig) — covering both pure `Virus` amplicons and human–viral hybrids. `True`/`False`. Only meaningful when the `GRCh38_viral` reference was used; `False` otherwise. |
| `BFBArchitect_min_score`       | Best BFBArchitect reconstruction score observed for the amplicon, or `NA` when no score is available. Lower scores indicate more parsimonious reconstructions. |
| `BFB_source`                   | Source of a positive BFB call: `AC`, `BFBArchitect`, or `AC\|BFBArchitect`; `NA` when no BFB call was made. |

When `--verbose_classification` is set, the profile additionally includes raw classification scores and extended BFBArchitect diagnostics such as passing-region count, multiplicities, attempted regions, and whole-graph/reverse-polarity status.

The `amplicon_decomposition_class` is an abstract label and can be one of six classes:

- `Cyclic`: This indicates the amplicon is bioinformatically cyclic (genome cycles) - and may be either an ecDNA or BFB (check `ecDNA+` and `BFB+` columns)
- `Complex non-cyclic`: (CNC) The amplicon contains a focal amplification with significant rearrangements (e.g. derived by chromothripsis), but does not contain genome cycles characteristic of ecDNA. However, this may class still contain a BFB (check `BFB+` column).
- `Linear`: A focal amplification with few to no significant rearrangments evident - frequently the exact mechanism is unclear. Label also includes low CN focal amplifications caused by tandem duplications.
- `No-FSCNA`: The AA amplicon is valid, but AC did not detect a focal significant copy-number amplification.
- `Invalid`: The AA amplicon failed AC validity checks, such as when cycles are removed by low-complexity filtering. LC filtering removes paths/cycles only when their LC-overlapping bp fraction or LC-overlapping discordant-breakend fraction exceeds configured thresholds. AC also marks high foldback-orientation artifact cases as `Invalid` and logs QC guidance for re-running AmpliconSuite-pipeline with stricter foldback pair support.
- `Virus`: If the GRCh38_viral reference was used, then this amplicon corresponds to an oncoviral genome.

#### ****`[prefix]_gene_list.tsv`****
Reports the genes present on amplicons with each classification, and which genomic feature (e.g. ecDNA_1, BFB_1, etc), it is located on, along with the copy number and which end(s) of the gene have been lost ("truncated"), will be one of `None`, `5p` (5-prime end), `3p` (3-prime end) or `5p_3p` if both. Genes are sourced from RefGene and most lncRNAs and micro-RNAs are excluded from the report.

 | Column name             | Contents                                                                                                                                 |
|-------------------------|------------------------------------------------------------------------------------------------------------------------------------------|
| `sample_name`           | Sample name prefix                                                                                                                       |
| `amplicon_number`       | AA amplicon index, e.g. `amplicon2`                                                                                                      |
| `feature`               | Which feature inside the amplicon the gene is present on. May be `unknown` if cannot be confidently assigned to a feature.               |
| `gene`                  | Gene name (RefGene)                                                                                                                      |
| `gene_cn`               | Maximum copy number of genomic segments (larger than 1kbp) overlapping the gene, as reported by AA                                       |
| `truncated`             | Which end(s) of the gene have been lost ("truncated"), will be one of `None`, `5p` (5-prime end), `3p` (3-prime end) or `5p_3p` if both  |
| `is_canonical_oncogene` | Reports if gene is present in [COSMIC](https://cancer.sanger.ac.uk/cosmic/curation), [ONGene](https://ongene.bioinfo-minzhao.org/).      |
| `ncbi_id`               | Reports the [NCBI Accession ID](https://www.ncbi.nlm.nih.gov/books/NBK470040/) of the gene                                               |

#### ****`[prefix]_lncRNA_list.tsv`****
This file has a highly similar structure to the `gene_list.tsv` file and is based on [GENCODE](https://www.gencodegenes.org/) lncRNA annotations. Note some genes overlap with GENCODE lncRNA and the genes. Those are primarily reported in the gene list file.

#### ****`[prefix]_feature_basic_properties.tsv`****
Reports a table of basic properties such as size of captured regions, median and max CN, and a flag field to report if the call is "borderline" (ecDNA with CN < 8, other classes with CN < 5).

#### ****`[prefix]_feature_complexity.tsv`****
Reports per-feature complexity scores as measured by the number of genomic segments and the diversity of copy number among all the amplicon decompositions performed by AA. For more information please see the Supplementary Information file of [this study](https://www.nature.com/articles/s41586-023-05937-5).

 | Column name                     | Contents   |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample_name`                   | Sample name prefix |
| `amplicon_number`               | AA amplicon index, e.g. `amplicon2` |
| `feature`                       | Which feature this score applies to, e.g. `ecDNA_1`. |
| `feature_complexity`            | The feature complexity score. This is the value reported as "Complexity score" in the results table. |
| `decomp_complexity`             | Amount of complexity or diversity captured in the AA decompositions overlapping this feature, without the segment-count term. |
| `nseg_complexity`               | Complexity contribution from the number of genomic segments overlapping this feature. |


#### ****`[output_prefix]_fan_calls.tsv`****
Reports the FAN (Focal amplification in neochromosome) classifier result for each AA amplicon. The file includes `sample_name`, `amplicon_number`, `fan_call` (`Positive`/`None detected`), `fan_probability`, and the model features used for the call.

**What FAN detects.** FAN identifies a focal-amplification phenotype that is structurally distinct from ecDNA. Where ecDNA is spatially focal, high copy-number, and circular, a FAN amplicon is broad and copy-number–heterogeneous, carries a high density of large or interchromosomal structural variants, and spreads across two or more chromosomes — a pattern that resolves into a chromosomal or neochromosomal (i.e. non-circular, integrated) structure rather than an extrachromosomal circle. A neochromosome is a large, structurally rearranged chromosome assembled *de novo* from amplified material, frequently drawing on multiple chromosomes of origin and carrying amplified oncogenes (Garsed et al., *[The architecture and evolution of cancer neochromosomes](https://doi.org/10.1016/j.ccr.2014.09.010)*, Cancer Cell, 2014). The FAN call is produced by a lightweight machine-learning model that scores each amplicon's AA graph from a compact set of features capturing SV burden, chromosomal spread, and amplified span; `fan_probability` is the model's score and `fan_call` applies the default decision threshold (0.5).

FAN is reported as an independent flag and is **not mutually exclusive** with `ecDNA+` or `BFB+` — a single amplicon can be both FAN+ and ecDNA+. In cell-line FISH validation, FAN calls corresponded predominantly to non-native HSRs (amplified material integrated at a non-native chromosomal location), consistent with the neochromosomal interpretation. No FAN call corresponded to a double minute (DM) alone: in the cases where a FAN-classified amplicon coincided with a FISH-observed DM, that amplicon was also called `ecDNA+`.

#### ****`[output_prefix]_feature_similarity_scores.tsv`****
Reports pairwise feature similarity scores for cross-sample feature pairs with overlapping genomic intervals. Same-sample pairs are not compared. The score table includes ecDNA, BFB, Linear, Complex-non-cyclic, and FAN features, and is written by default for each AC run. These rows are also the audit trail used by `--filter_similar`; the table is written regardless of whether filtering is enabled.

#### ****`[output_prefix]_ecDNA_counts.tsv`****
This two-column file reports the `sample_name` and the number of ecDNA identified in the samples.

#### ****`[output_prefix]_ecDNA_context_calls.tsv`****
This two column file reports the ecDNA feature name (sample_amplicon_ecDNA_number), and a classification of the ecDNA focal amplification genome context.
The possibilities for ecDNA context classification are

| ecDNA context class                         | Description                                                                                                                                                                                         |
|-------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Simple circular simple background   | A simple ecDNA cycle with minimal rearrangements in the surrounding genome. Likely not derived from chromothripsis.                                                                              |
| Simple circular complex background  | A simple ecDNA cycle however there are genomic rearrangements in the vicinity outside the ecDNA region.                                                                                          |
| BFB-like | ecDNA possibly derived from a BFB.                                                                                                                                                               |
| Two-foldback | ecDNA being flanked by two foldback-like SVs. Likely not derived from chromothripsis, but possibly from [ODIRA](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002016). |
| Heavily rearranged unichromosomal | ecDNA from a heavily rearranged genome on one chromosome. Possibly due to chromothripsis. |
| Heavily rearranged multichromosomal | ecDNA from a heavily rearranged genome involving multiple chromosomes. Possibly due to chromothripsis and chromoplexy. |
| Unknown | Does not match any of the classes above. |

We suggest that `Simple circular simple background` and `Two-foldback` are most closely associated with simple excisional ecDNA. `Simple circular complex background`, `BFB-like`, `Heavily rearranged uni/multichromosomal` are most closely associated with chromothripsis.

#### Amplicon bed files, annotated cycles, and SV summaries
Additionally, there are three directories  created by `amplicon_classifier.py`. They are
- `[prefix]_classification_bed_files/`, which contains bed files of the regions classified into each feature. May contain bed files marked `unknown` if the region could not be confidently assigned.
  - FAN-positive amplicons are reported as a numbered feature, e.g. `sample_amplicon1_FAN_1_intervals.bed`.
  - The bed files report genomic intervals using a [0-based, half-open counting system](https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/). This is the same system used by the UCSC genome browser.
  - By contrast, AmpliconArchitect's graph and cycles files report genomic coordinates using a 0-based, fully closed counting system. This means that intervals reported by AC will contain one additional base on the second coordinate, which is not part of the amplicon (half-open).
  - Intervals reported in these bed files do not represent the structure of ecDNA, and may face limitations related to missing SVs or inexactly refined amplicon endpoints (a limitation of short-reads).

- `[prefix]_SV_summaries/`, which contains tab-separated files summarizing the SVs detected by AA and what features the overlap in the amplicon.
- `[prefix]_annotated_cycles_files/`, which contains AA cycles files with additional annotations about length of discovered paths/cycles and their classification status. Using these annotated cycles files is preferred over the unannotated cycles file produced by AA. These cycles are filtered to remove paths/cycles with excessive low-complexity content, patch reference genome issues, and filter duplicate cycle entries erroneously output by AA (uncommon).

#### Results table
By default (in `--AA_results`/`--input` mode), AC creates `[prefix]_result_table.tsv` and `[prefix]_result_data.json`; pass `--no_results_table` to skip this. The table has one row per reported feature, including ecDNA, BFB, decomposition features, and FAN features. FAN rows use the numbered feature ID `FAN_1`; the `FAN probability` column reports the model probability for each amplicon, and the `Contains viral` column carries the amplicon's `contains_viral` status.


### 4. Command-Line Options

If running AC only on a single AA amplicon, use arguments:
- `-c/--cycles`: AA cycles file
- `-g/--graph`: AA graph file

Else if running on multiple amplicons, use argument
- `--AA_results`: Path to a directory containing one or more AA results. AC will search this directory recursively and index all the AA results it finds for classification.


| Column name                                     | Contents                                                                                                                                                                                                             |
|-------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| 
| `--ref [hg19, GRCh37, GRCh38, GRCh38_viral, mm10, or GRCm38]` | (Required) Choose reference genome version used in generating AA output.                                                                                                                                |
| `-v/--version`                                  | Print version and exit.                                                                                                                                                                                              |
| `-o`                                            | Output filename prefix. Can include a custom dir path. Defaults to prefix of `-i` or `-c`.                                                                                                                           |
| `--add_chr_tag`                                 | Adds back missing "chr" prefix to chromosome names. If you have a mix of hg19 and GRCh37 amplicons, set `--ref hg19` and `--add_chr_tag` to classify them all together.                                              |
| `--min_flow`                                    | Minumum cycle CN flow to consider among decomposed paths (default=1).                                                                                                                                                | 
| `--min_size`                                    | Minimum cycle size (in bp) to consider as valid amplicon (default=5000).                                                                                                                                             |
| `--verbose_classification`                      | Add raw classification scores and extended BFBArchitect diagnostics to the `amplicon_classification_profiles.tsv` file. `BFBArchitect_min_score` and `BFB_source` are included by default. Useful for debugging.              |
| `--report_complexity`                           | [Deprecated — on by default] Compute an amplicon entropy/complexity measure for each amplicon (written to `_feature_complexity.tsv`).                                                                                |
| `--force`                                       | Disable No-FSCNA/Invalid class, if possible. Use only when extremely large CN seeds were used in AA amplicon generation (>10 Mbp intervals) or if debugging.                                                           |
| `--decomposition_strictness`                    | Value between 0 and 1 reflecting how strictly to filter low CN decompositions (default = 0.1). Higher values filter more of the low-weight decompositions.                                                           |
| `--exclude_bed`                                 | Provide a bed file of regions to ignore during classification. Useful for separating linked amplicons or augmenting existing low-complexity annotations.                                                             |
| `--no_LC_filter`                                | Set this to turn off filtering low-complexity & poor mappability genome region paths & cycles based on the regions in the AA data repo.                                                                              |
| `--filter_similar`                              | Permits filtering of false-positive amps arising in multiple independent samples based on similarity calculation. Only use if all samples are of independent origins (not replicates and not multi-region biopsies). |
| `--filter_pval`                                 | P-value cutoff used when `--filter_similar` is set. Default is `0.05/(n_amps-1)`.                                                                                                                                    |
| `--config`                                      | Path to a custom parameter configuration file. If not specified, uses the default config in the `ampclasslib` directory.                                                                                            |
| `--no_results_table`                            | Skip creation of the summary results table. The `_result_table.tsv`/`_result_data.json` table is generated by default (when running in `--AA_results`/`--input` mode); use this flag to turn it off.                   |
| `--no_bfbarchitect`                             | Disable BFBArchitect integration. By default, AmpliconClassifier uses BFBArchitect when it is installed and warns when it is unavailable.                                                                               |
| `--bfb_threads`                                 | Number of BFBArchitect ILP solver threads per amplicon. For multi-amplicon runs, `2` is generally the best efficiency setting. When using `--jobs`, approximate total BFBArchitect solver threads are `--jobs * --bfb_threads`. |
| `--jobs`                                        | Number of amplicons to classify in parallel. Default is 1. Output rows are collected in input order.                                                                                                                   |
| `-i/--input`                                    | If you have already run `make_input.sh`, you can give the resulting .input file instead of setting `--AA_results`                                                                                                    | 

### 5. Other utilities:

### Feature Similarity
One may wish to compare two overlapping focal amplifications and quantify their similarity - particularly when they are derived from multi-region 
or longitudinal sampling. We provide a script which ***a)*** identifies overlap between pairs of amplicons (using the same input file as `amplicon_classifier.py`), 
***b)*** computes measurements of the similarity of the two overlapping amplicons based on shared breakpoints and shared genomic content - 
using both a Jaccard index approach and also our own *Symmetric Similarity Score* and *Asymmetric Similarity Score* approaches, and ***c)*** compares the scores against
the similarity scores for overlapping amplicons derived from unrelated origins (data derived from Turner et al. _Nature_ 2017 and deCarvalho et al. _Nature Genetics_ 2018, Bergstrom et al. _Nature_ 2020 and Paulson et al. _Nature Communications_ 2022).
The output file `*_similarity_scores.tsv` reports the following columns:
- Amplicon 1 ID & Amplicon 2 ID
- Symmetric Similarity Score (a combination of GenomicSegment and Breakpoint scores)
- Percentile and P-value of the Sym. Score in the background unrelated overlapping amplicon distribution. P-value based on beta distribution fit to similarity scores.
- `GenomicSegmentScore1` & `GenomicSegmentScore2` based on the directional similarity of genomic segment overlap (Amp1 and Amp2)/(Amp1) or (Amp1 and Amp2)/(Amp2), respectively.
- `BreakpointScore1` & `BreakpointScore2` based on the directional similarity of breakpoint matching (Amp1 and Amp2)/(Amp1) or (Amp1 and Amp2)/(Amp2), respectively.
- `JaccardGenomicSegment`, based on overlap of genomic segments (based on overlap of genomic coordinates)
- `JaccardBreakpoint`, based on overlap from matching of breakpoints.

The primary difference between `amplicon_similarity.py` and `feature_similarity.py` is that the latter only considers regions annotated as specific classification features (e.g. ecDNA, BFB) when computing similarity scores.

**Example command for `feature_similarity.py`**
```shell
$AC_SRC/feature_similarity.py --ref hg38 -f [sample]_features_to_graph.txt --required_classifications [default "any", but can be "ecDNA", "BFB", "CNC", "linear" or any combination of those] 
```
Where "[sample]_features_to_graph.txt" is one of the output files generated by `amplicon_classifier.py` and is a two-column text file with the feature ID (e.g. `sample_amplicon1_ecDNA1`) as the first column and the path of the associated AA graph file in the second column. `--include_path_in_feature_name` is useful for comparing two runs with the same sample names 
(after combining both `_features_to_graph.txt` files from each run into one file).

**How do I use this to perform similarity score-based filtering on my samples?**
To remove potential false-positive focal amplification calls from a collection of samples, AC uses the similarity scores and p-value reported in `[output_prefix]_feature_similarity_scores.tsv`. The motivating idea is that in unrelated samples, precisely conserved CN boundaries and SVs should be exceptionally rare unless they are derived from issues with the reference genome. If all samples in the collection are from unrelated sources (no replicates, no multi-region or longitudinal samples, etc.) users can simply run AC on the batch of samples setting the `--filter_similar` flag. AC never compares features from the same sample during this filtering step.

When `--filter_similar` is set, AC reuses the same feature similarity rows that it writes to `[output_prefix]_feature_similarity_scores.tsv`. ecDNA, BFB, Linear, and Complex-non-cyclic hits are filtered at feature scope. FAN hits are filtered at amplicon scope: if a FAN feature is significantly similar to a feature from another sample, AC removes the whole amplicon from reported feature outputs and downgrades the amplicon to `No-FSCNA`. The FAN call for that amplicon is also cleared to `None detected`. The triggering similarity rows remain in the similarity score table for review.

If the collection contains a mixture of samples from related and unrelated sources, then some samples will likely have focal amplifications having a high degree of similarity due to being from related origins. In this case, users can run the `feature_similarity.py` script described above on the samples to produce a table of similarity scores and p-values. Users can then filter highly similar focal amplifications from unrelated samples using the p-value in the table. Since this involves multiple hypothesis testing, to control false positive rate, users can mirror what is done internally by AC and apply a slightly modified Bonferroni correction of `alpha/(n-1)` where alpha by default is 0.05 and `n` is the number of samples in the collection having a focal amplification.

### SV Support Checker


Counts read support for calls from a AA SV summary TSV table by checking for supporting read evidence in a BAM file.

#### What it does
- Takes SV calls (typically from AmpliconArchitect) and searches a BAM file for discordant read pairs and split reads that support each SV
- Matches reads to SVs based on breakpoint positions (within configurable tolerance)
- Considers both position AND orientation to ensure reads truly support the SV topology
- Outputs the original table with three additional columns showing support counts from the BAM

#### Usage
```bash
check_SV_support.py \
  --input_table sv_calls.tsv \
  --bam sample.bam \
  --ref GRCh38
```

#### Common Use Cases
- Confirm raw support for SV calls in a sample using the BAM
- Check if SVs detected in sample are supported by reads in a different sample (e.g. find support for rare ecDNA in a related sample)

#### Output
The script outputs the original TSV with three additional columns:
- `bam_discordant_pairs` - Number of discordant read pairs supporting the SV
- `bam_split_reads` - Number of split reads supporting the SV  
- `bam_total_support` - Total support (sum of above)

Output file is written to the current working directory with the format:
```
{input_name}_annotated_{bam_name}.tsv
```

#### Options
- `--input_table, -i` (required) - Input TSV file containing SV calls
- `--bam, -b` (required) - BAM file to validate SVs against
- `--ref, -r` (required) - Reference genome (hg19, GRCh37, GRCh38, hg38, mm10)
- `--output, -o` (optional) - Custom output file path
- `--tolerance, -t` (default: 150bp) - Breakpoint position tolerance, can be increased for less strict matching
- `--window, -w` (default: 500bp) - Search window around each breakpoint
- `--min_mapq` (default: 5) - Minimum mapping quality
- `--min_sv_size` (default: 0) - Minimum SV size for intrachromosomal events
- `--verbose, -v` - Show detailed orientation matching for debugging

#### Example
```bash
# Basic usage
check_SV_support.py -i COLO320_amplicon3_SV_summary.tsv \
                     --bam COLO320_metastasis.bam \
                     --ref GRCh38

# With custom tolerance and verbose output
check_SV_support.py -i sv_calls.tsv \
                     --bam sample.bam \
                     --ref hg38 \
                     --tolerance 200 \
                     --verbose
```

#### Notes
- Requires AA_DATA_REPO environment variable to be set (for low complexity regions)
- Input table should have standard SV format with columns: chrom1, pos1, chrom2, pos2, sv_type, orientation, etc.
- Only counts reads where both position AND orientation match the SV
- If SV orientation is not available, counts based on position matching only


### Info about accessing the legacy version of AmpliconClassifier:

For the legacy version used in [Kim et al., *Nature Genetics*, 2020](https://www.nature.com/articles/s41588-020-0678-2) please see the scripts and README in the "legacy_natgen_2020" folder of this repo.
The legacy version is only recommended for reproducing the paper results, and not for state-of-the-art amplicon classification. The legacy version was developed by Nam Nguyen, Jens Luebeck, and Hoon Kim.
The current version was developed by Jens Luebeck and Bhargavi Dameracharla.
