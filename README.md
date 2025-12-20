# AmpliconClassifier
![GitHub release (latest by date)](https://img.shields.io/github/v/release/jluebeck/AmpliconClassifier?display_name=release)
![GitHub](https://img.shields.io/github/license/jluebeck/AmpliconClassifier)

### Classify [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) outputs to predict types of focal amplifications present.

This tool classifies the outputs of [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconSuite-pipeline).

If using AmpliconClassifier, please cite the following publication which describes the AmpliconClassifier methodology in the [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05937-5/MediaObjects/41586_2023_5937_MOESM1_ESM.pdf) section:

Luebeck J, et al. [Extrachromosomal DNA in the cancerous transformation of Barrett’s oesophagus](https://www.nature.com/articles/s41586-023-05937-5). Nature. 2023. PMID: 37046089

See the section at the end of the README for information about the legacy version.<br />
<br />

### 1. Installation

AmpliconClassifier is included with [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline), but for re-classification, you may wish to keep a standalone installation of this module, using the instructions below.

#### Requirements
- Either python 2.7+ or python 3.0+
- `intervaltree`, `scipy`, `pandas`, `pysam` python libraries:  
>`conda install intervaltree scipy pandas pysam`  # or `pip install intervaltree scipy pandas pysam`
- `$AA_DATA_REPO` environment variable and data repo files. See instructions [here](https://github.com/AmpliconSuite/AmpliconArchitect#setting-up-the-aa-data-repo). 
- (optional) `matplotlib`: `conda install matplotlib-base` or `pip install matplotlib`

Mac users will need to perform one additional installation step:
```bash
brew install coreutils
```

#### Setup:
```bash
git clone https://github.com/jluebeck/AmpliconClassifier.git
cd AmpliconClassifier
echo export AC_SRC=$PWD >> ~/.bashrc
source ~/.bashrc
```


### 2. Usage

`amplicon_classifier.py` takes a collection of (or single) AA graph files and corresponding AA cycles file as inputs.

**Most common - classifying multiple amplicons:**
You can provide the directory containing multiple AA amplicons or multiple uniquely named samples
>`python amplicon_classifier.py --ref GRCh38 --AA_results /path/to/AA/output/directories/ > classifier_stdout.log`

AC will crawl the given location and find all relevant AA files and perform classification on them.

**To classify a single amplicon**:

>`python amplicon_classifier.py --ref GRCh38 --cycles sample_amplicon1_cycles.txt --graph sample_amplicon1_graph.txt > classifier_stdout.log`


**Less common - separate usage of `make_input.sh`:**

Alternatively, you can use the `make_input.sh` script to gather the necessary input files outside of AC:

`make_input.sh` takes a path and an output prefix. e.g:

>`make_input.sh /path/to/AA/output/directories/ example_collection` 

This would create a file called `example_collection.input` which can be given as the `--input` argument for AC.


#### Combining classification results from GRCh37 and hg19:

If combining data from both GRCh37 and hg19 in the same classification run, you can set the flag `--add_chr_tag` to add the "chr" prefix to each chromosome name and effectively unify everything as hg19-based.

### 3. Outputs

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

#### ****`[prefix]_feature_entropy.tsv`****
Reports amplicon complexity scores as measured by the number of genomic segments and the diversity of copy number among all the amplicon decompositions performed by AA. For more information please see the Supplementary Information file of [this study](https://www.nature.com/articles/s41586-023-05937-5).

 | Column name                     | Contents   |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample_name`                   | Sample name prefix |
| `amplicon_number`               | AA amplicon index, e.g. `amplicon2` |
| `feature`                       | Which feature inside the amplicon the gene is present on. May be `unknown` if cannot be confidently assigned to a feature. |
| `total_feature_entropy`         | This is the amplicon complexity score. |
| `decomp_entropy`                | Amount of entropy or diversity captured in the AA decompositions overlapping this feature. |
| `Amp_nseg_entropy`              | Amount of entropy or diversity captured by the number of genomic segments overlapping this feature. |


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
  - The bed files report genomic intervals using a [0-based, half-open counting system](https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/). This is the same system used by the UCSC genome browser.
  - By contrast, AmpliconArchitect's graph and cycles files report genomic coordinates using a 0-based, fully closed counting system. This means that intervals reported by AC will contain one additional base on the second coordinate, which is not part of the amplicon (half-open).
  - Intervals reported in these bed files do not represent the structure of ecDNA, and may face limitations related to missing SVs or inexactly refined amplicon endpoints (a limitation of short-reads).

- `[prefix]_SV_summaries/`, which contains tab-separated files summarizing the SVs detected by AA and what features the overlap in the amplicon.
- `[prefix]_annotated_cycles_files/`, which contains AA cycles files with additional annotations about length of discovered paths/cycles and their classification status. Using these annotated cycles files is preferred over the unannotated cycles file produced by AA. These cycles are filtered to remove cycles overlapping low-complexity regions of the genome, patches reference genome issues, and filters duplicate cycle entries erroneously output by AA (uncommon).


### 4. Command-Line Options

If running AC only on a single AA amplicon, use arguments:
- `-c/--cycles`: AA cycles file
- `-g/--graph`: AA graph file

Else if running on multiple amplicons, use argument
- `--AA_results`: Path to a directory containing one or more AA results. AC will search this directory recursively and index all the AA results it finds for classification.


| Column name                                     | Contents                                                                                                                                                                                                             |
|-------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| 
| `--ref [hg19, GRCh37, GRCh38, mm10, or GRCm38]` | (Required) Choose reference genome version used in generating AA output.                                                                                                                                             |
| `-v/--version`                                  | Print version and exit.                                                                                                                                                                                              |
| `-o`                                            | Output filename prefix. Can include a custom dir path. Defaults to prefix of `-i` or `-c`.                                                                                                                           |
| `--add_chr_tag`                                 | Adds back missing "chr" prefix to chromosome names. If you have a mix of hg19 and GRCh37 amplicons, set `--ref hg19` and `--add_chr_tag` to classify them all together.                                              |
| `--min_flow`                                    | Minumum cycle CN flow to consider among decomposed paths (default=1).                                                                                                                                                | 
| `--min_size`                                    | Minimum cycle size (in bp) to consider as valid amplicon (default=5000).                                                                                                                                             |
| `--verbose_classification`                      | Output verbose information in the `amplicon_classification_profiles.tsv` file, and create `edge_classification_profiles.tsv`. Useful for debugging.                                                                  |
| `--force`                                       | Disable No amp/Invalid class, if possible. Use only when extremely large CN seeds were used in AA amplicon generation (>10 Mbp intervals) or if debugging.                                                           |
| `--plotstyle [noplot, individual]`              | \[experimental] Produce a radar-style plot of classification strengths. Default `noplot`.                                                                                                                            |
| `--decomposition_strictness`                    | Value between 0 and 1 reflecting how strictly to filter low CN decompositions (default = 0.1). Higher values filter more of the low-weight decompositions.                                                           |
| `--exclude_bed`                                 | Provide a bed file of regions to ignore during classification. Useful for separating linked amplicons or augmenting existing low-complexity annotations.                                                             |
| `--no_LC_filter`                                | Set this to turn off filtering low-complexity & poor mappability genome region paths & cycles based on the regions in the AA data repo.                                                                              |
| `--filter_similar`                              | Permits filtering of false-positive amps arising in multiple independent samples based on similarity calculation. Only use if all samples are of independent origins (not replicates and not multi-region biopsies). |
| `--make_results_table`                          | Creates summary results table (_results_table.tsv) after classification completes.                                                                                                                                   |
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
Where "[sample]_features_to_graph.txt" is one of the output files generated by `amplicon_classifier.py`. `--include_path_in_features` is useful for comparing two runs with the same sample names 
(after combining both `_features_to_graph.txt` files from each run into one file).

**How do I use this to perform similarity score-based filtering on my samples?**
To remove potential false-positive focal amplificaiton calls from a collection of samples, we can use the similarity scores and p-value reported by AC. The motivating idea is that in unrelated samples, precisely conserved CN boundaries and SVs should be exceptionally rare unless they are derived from issues with the reference genome. If all samples in the collection are from unrelated sources (no replicates, no multi-region or longitudinal samples, etc.) users can simply run AC on the batch of samples setting the `--filter_similar` flag. However, if the collection contains a mixture of samples from related and unrelated sources, then some samples will likely have focal amplifications having a high degree of similarity due to being from related origins. In this case, users can run the `feature_similarity.py` script described above on the samples to produce a table of similarity scores and p-values. Users can then filter highly similar focal amplifications from unrelated samples using the p-value in the table. Since this involves multiple hypothesis testing, to control false positive rate, users can mirror what is done by internally by AC and apply a slightly modified Bonferroni correction of `alpha/(n-1)` where alpha by default is 0.05 and `n` is the number of samples in the collection having a focal amplification. 

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
