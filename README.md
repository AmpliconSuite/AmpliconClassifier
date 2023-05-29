# AmpliconClassifier
![GitHub release (latest by date)](https://img.shields.io/github/v/release/jluebeck/AmpliconClassifier?display_name=release)
![GitHub commits since latest release (by date)](https://img.shields.io/github/commits-since/jluebeck/AmpliconClassifier/v0.5.1/main)
![GitHub](https://img.shields.io/github/license/jluebeck/AmpliconClassifier)

### Classify [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) outputs to detect types of focal amplifications present.
### Current version: 0.5.1

This tool classifies the outputs of AmpliconArchitect, [which is available here](https://github.com/jluebeck/AmpliconArchitect).

If using AmpliconClassifier (current version, not legacy version), please cite:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Luebeck et al., [Extrachromosomal DNA in the cancerous transformation of Barrett's esophagus](https://www.nature.com/articles/s41586-023-05937-5).
*Nature*. 2023.

If referring to the legacy version, see the section at the end of the README.
<br />
<br />

***Please note that this software is actively being developed. Stable versions are released on the main branch.***

### 1. Prerequisites: 
- Either python2 or python3.
- `intervaltree` (python library):  `pip install intervaltree`
- `scipy` (python library): `pip install scipy`
- `$AA_DATA_REPO` environment variable. For instructions see [AmpliconArchitect installation](https://github.com/jluebeck/AmpliconArchitect#data-repositories). 
- `matplotlib` (python library, optional): `pip install matplotlib`

#### Installation:
```bash
git clone https://github.com/jluebeck/AmpliconClassifier.git
cd AmpliconClassifier
echo export AC_SRC=$PWD >> ~/.bashrc
source ~/.bashrc
```

Mac users will need to perform one additional installation step:
```bash
brew install coreutils
```

### 2. Usage:

`amplicon_classifier.py` takes a collection of (or single) AA graph files and corresponding AA cycles file as inputs.

To classify a single amplicon,

`python amplicon_classifier.py --ref [hg19, GRCh37, or GRCh38] --cycles sample_amplicon1_cycles.txt --graph sample_amplicon1_graph.txt > classifier_stdout.log`

**More commonly, if classifying multiple amplicons, you can use the `make_input.sh` script to gather the necessary input files automatically.**

`make_input.sh` takes a path and an output prefix. e.g:

`make_input.sh /path/to/AA/output/directories/ example_collection` 

This will create a file called `example_collection.input` which can be given as the `--input` argument for AC.

To subsequently generate classifications for a list of amplicons:

`python amplicon_classifier.py --ref [hg19, GRCh37, or GRCh38] --input [file with list of your amplicons] > classifier_stdout.log`

and it will search for the cycles and graph files in that directory, and pair the locations into a text file compatible with the `--input` argument.

There is also an experimental option you can set to visualize the strength of each amplicon class assigned to an amplicon, which can be turned on by setting `--plotStyle individual`.

If combining data from both GRCh37 and hg19 in the same classification run, you can set the flag `--add_chr_tag` to add the "chr" prefix to each chromosome name and effectively unify everything as hg19-based.

### 3. Output:

#### ****`[prefix]_amplicon_classification_profiles.tsv`**** 

Contains an abstract classification of the amplicon, and also indicates in separate columns "BFB+" and "ecDNA+" status.
Note that amplicons receiving a "Cyclic" classification may be ecDNA+, BFB+ or both.

| Column name                    | Contents                                                                                                                                                                                |
|--------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample_name`                  | Sample name prefix                                                                                                                                                                      |
| `amplicon_number`              | AA amplicon index, e.g. `amplicon2`                                                                                                                                        |
| `amplicon_decomposition_class` | Abstract description of the AA amplicon type. Note that `Cyclic` can refer to either BFB or ecDNA. Please see the following columns for that distinction.                                                                                    |
| `ecDNA+`                       | Prediction about whether the AA amplicon contains ecDNA. Note, an AA amplicon may contain regions surrounding the ecDNA, or multiple linked ecDNA. Either `Positive` or `None detected` |
| `BFB+`                         | Prediction about whether the AA amplicon is the result of a BFB. Either `Positive` or `None detected`                                                                                   |
| `ecDNA_amplicons`              | Predicted number of distinct (non-overlapping) ecDNA which are represented in a single AA amplicon. This estimate is highly experimental.                                               |

Because an ecDNA may overlap with a BFB, they are reported separately.

#### ****`[prefix]_gene_list.tsv`****
Reports the genes present on amplicons with each classification, and which genomic feature (e.g. ecDNA_1, BFB_1, etc), it is located on, along with the copy number and which end(s) of the gene have been lost ("truncated"), will be one of `None`, `5p` (5-prime end), `3p` (3-prime end) or `5p_3p` if both. Genes are sourced from RefGene and most lncRNAs and micro-RNAs are excluded from the report.

 | Column name                    | Contents                                                                                                                                                                                                                                                                   |
|--------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample_name`                  | Sample name prefix                                                                                                                                                                                                                                                         |
| `amplicon_number`              | AA amplicon index, e.g. `amplicon2`                                                                                                                                                                                                                                        |
| `feature` | Which feature inside the amplicon the gene is present on. May be `unknown` if cannot be confidently assigned to a feature.                                                                                                                                                 |
| `gene`                       | Gene name (RefGene)                                                                                                                                                                                                                                                        |
| `gene_cn`                         | Maximum copy number of genomic segments (larger than 1kbp) overlapping the gene, as reported by AA                                                                                                                                                                         |
| `truncated`              | Which end(s) of the gene have been lost ("truncated"), will be one of `None`, `5p` (5-prime end), `3p` (3-prime end) or `5p_3p` if both                                                                                                                                    |
| `is_canonical_oncogene` | Reports if gene is present in [COSMIC](https://cancer.sanger.ac.uk/cosmic/curation), [ONGene](https://ongene.bioinfo-minzhao.org/). |

#### ****`[prefix]_feature_basic_properties.tsv`****
Reports a table of basic properties such as size of captured regions, median and max CN, and a flag field to report if the call is "borderline" (ecDNA with CN < 8, other classes with CN < 5).

#### ****`[prefix]_feature_entropy.tsv`****
Reports amplicon complexity scores as measured by the number of genomic segments and the diversity of copy number among all the amplicon decompositions performed by AA. For more information please see the Supplementary Information file of [this study](https://www.nature.com/articles/s41586-023-05937-5).

 | Column name                    | Contents   |
|--------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample_name`                  | Sample name prefix |
| `amplicon_number`              | AA amplicon index, e.g. `amplicon2` |
| `feature` | Which feature inside the amplicon the gene is present on. May be `unknown` if cannot be confidently assigned to a feature. |
| `total_feature_entropy`                       | This is the amplicon complexity score. |
| `decomp_entropy`                         | Amount of entropy or diversity captured in the AA decompositions overlapping this feature. |
| `Amp_nseg_entropy`              | Amount of entropy or diversity captured by the number of genomic segments overlapping this feature. |


#### ****`[output_prefix]_ecDNA_counts.tsv`****
This two-column file reports the `sample_name` and the number of ecDNA identified in the sample across all amplicons from the sample.

#### Amplicon bed files, annotated cycles, and SV summaries
Additionally, there are three directories  created by `amplicon_classifier.py`. They are
- `[prefix]_classification_bed_files/`, which contains bed files of the regions classified into each feature. May contain bed files marked `unknown` if the region could not be confidently assigned.
  - The bed files report genomic intervals using a [0-based, half-open counting system](https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/). This is the same system used by the UCSC genome browser.
  - By contrast, AmpliconArchitect's graph and cycles files report genomic coordinates using a 0-based, fully closed counting system. This means that intervals reported by AC will contain one additional base on the second coordinate, which is not part of the amplicon (half-open).
  - Intervals reported in these bed files do not represent the structure of ecDNA, and may face limitations related to missing SVs or inexactly refined amplicon endpoints (a limitation of short-reads).

- `[prefix]_SV_summaries/`, which contains tab-separated files summarizing the SVs detected by AA and what features the overlap in the amplicon.
- `[prefix]_annotated_cycles_files/`, which contains AA cycles files with additional annotations about length of discovered paths/cycles and their classification status. Using these annotated cycles files is preferred over the unannotated cycles file produced by AA.


### 4. Description of command line arguments:

Running on a single AA amplicon:
- `-c/--cycles`: AA cycles file
- `-g/--graph`: AA graph file

OR running on multiple amplicons
- `-i/--input`: Tab-separated file containing one or more amplicons formatted as

`sample_name_amplicon1   /path/to/sample_name_amplicon1_cycles.txt   /path/to/sample_name_amplicon1_graph.txt`

Other arguments
- `--ref [hg19, GRCh37, GRCh38, mm10, or GRCm38]`: (Required) Choose reference genome version used in generating AA output.
- `-v/--version`: Print version and exit.
- `-o`: Output filename prefix. Default is prefix of `-i` or `-c`.
- `--add_chr_tag`: If you have a mix of hg19 and GRCh37 amplicons, you can set `--ref hg19` and `--add_chr_tag` to classify them all together.
- `--min_flow`: Minumum cycle CN flow to consider among decomposed paths (default=1).
- `--min_size`: Minimum cycle size (in bp) to consider as valid amplicon (default=5000).
- `--report_comlexity`: Report a measurement of the amplicon's 'complexity' score, which represents a measurement of the complexity of the AA breakpoint graph decomposition.
- `--verbose_classification`: Output verbose information in the `amplicon_classification_profiles.tsv` file, and create `edge_classification_profiles.tsv`. Useful for debugging.
- `--force`: Disable No amp/Invalid class, if possible. Use only when extremely large CN seeds were used in AA amplicon generation (>10 Mbp intervals).
- `--plotstyle [noplot, individual]`: Produce a radar-style plot of classification strengths. Default `noplot`.
- `--decomposition_strictness`: Value between 0 and 1 reflecting how strictly to filter low CN decompositions (default = 0.1). Higher values filter more of the low-weight decompositions.
- `--exclude_bed`: Provide a bed file of regions to ignore during classification. Useful for separating linked amplicons or augmenting low-complexity annotations.
- `--no_LC_filter`: Set this to turn off filtering low-complexity & poor mappability genome region paths & cycles.
- `--no_LC_filter`: Set this to turn off filtering low-complexity & poor mappability genome region paths & cycles.
- `"--filter_similar`: Only use if all samples are of independent origins (not replicates and not multi-region biopsies). Permits filtering of false-positive amps arising in multiple independent samples based on similarity calculation.

### 5. Other utilities:

#### Amplicon Similarity and Feature Similarity
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

**Example command for `amplicon_similarity.py`**

```shell
$AC_SRC/amplicon_similarity.py --ref hg19 --add_chr_tag -i examples.input -o examples [--subset_bed some_intervals.bed] [--classification_file examples_amplicon_classification_profiles.tsv] [--required_classifications ecDNA BFB CNC]
```
Where `"examples.input"` is the tab-separate file generated by `make_input.sh /path/to/AA/outputs/ examples` and `--subset_bed` is an optional bed file restricting the similarity
calculations to one or more regions of the genome. Can provide classification file and list of classes to only perform similarity classifications on some amplicon classes.

**Example command for `feature_similarity.py`**
```shell
$AC_SRC/feature_similarity.py --ref hg38 -f [sample]_features_to_graph.txt --required_classifications [default "any", but can be "ecDNA", "BFB", "CNC", "linear" or any combination of those] 
```
Where "[sample]_features_to_graph.txt" is one of the output files generated by `amplicon_classifier.py`. `--include_path_in_features` is useful for comparing two runs with the same sample names 
(after combining both `_features_to_graph.txt` files from each run into one file).

**Info about accessing the legacy version of AmpliconClassifier:**

For the legacy version used in [Kim et al., *Nature Genetics*, 2020](https://www.nature.com/articles/s41588-020-0678-2) please see the scripts and README in the "legacy_natgen_2020" folder of this repo.
The legacy version is only recommended for reproducing the paper results, and not for state-of-the-art amplicon classification. The legacy version was developed by Nam Nguyen, Jens Luebeck, and Hoon Kim.
The current version is developed and maintained by Jens Luebeck.
