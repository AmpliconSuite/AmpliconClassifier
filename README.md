## AmpliconClassifier

#### Classify AmpliconArchitect output to detect types of focal amplifications present.

**Info about accessing the legacy version:** For the legacy version used in Kim et al., *Nature Genetics*, 2020 please see the scripts and README in the "legacy_natgen_2020" folder.
The legacy version is only recommended for reproducing the paper results, and not for state-of-the-art amplicon classification. The legacy version was developed by Nam Nguyen, Jens Luebeck, and Hoon Kim.
The current version is developed and maintained by Jens Luebeck.

If using AmpliconClassifier, please cite:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Kim H, Nguyen N, et al. [“Extrachromosomal DNA is associated with oncogene amplification and poor outcome across multiple cancers.”](https://www.nature.com/articles/s41588-020-0678-2)
*Nature Genetics*. 2020.

**Current version: 0.2.5 (stable)**

***Please note that this software is actively being developed. Stable versions are released on the main branch.***

**1. Prerequisites**: 
- Supports both python2 and python3.
- intervaltree (python library):  `pip install intervaltree`
- `$AA_DATA_REPO` environment variable. For instructions see [AmpliconArchitect installation](https://github.com/jluebeck/AmpliconArchitect#data-repositories). 
- matplotlib (python library, optional): `pip install matplotlib`

**2. Usage**:

`amplicon_classifier.py` takes an AA graph file and an AA cycles file as input.

To classify a single amplicon,

`python amplicon_classifier.py --ref [hg19, GRCh37, or GRCh38] --cycles [/path/to/amplicon_cycles.txt] --graph [/path/to/amplicon_graph.txt] > classifier_stdout.log`

Alternatively you can generate classifications for a list of amplicons:

`python amplicon_classifier.py --ref [hg19, GRCh37, or GRCh38] --input [file with list of your amplicons] > classifier_stdout.log`

If passing a list of amplicons, the `--input` argument must be formatted as follows, with one amplicon per line:

`sample_name_amplicon1   /path/to/sample_name_amplicon1_cycles.txt   /path/to/sample_name_amplicon1_graph.txt`

There is also an experimental option you can set to visualize the strength of each amplicon class assigned to an amplicon, which can be turned on by setting `--plotStyle individual`.

If you have data from both GRCh37 and hg19 that you are combining, you can set the flag `--add_chr_tag` to add the "chr" prefix to each chromosome name and effectively unify everything as hg19-based.

**3. Output**:

The most important file will be the file `[output_prefix]_amplicon_classification_profiles.tsv`. 
This contains an abstract classification of the amplicon, and also indicates in separate columns "BFB+" and "ecDNA+" status.
Note that amplicons receiving a "Cyclic" classification may be ecDNA+, BFB+ or both.

**4. Description of command line arguments**:

Running on a single AA amplicon:
- `-c/--cycles`: AA cycles file
- `-g/--graph`: AA graph file

OR running on multiple amplicons
- `-i/--input`: Tab-separated file containing one or more amplicons formatted as

`sample_name_amplicon1   /path/to/sample_name_amplicon1_cycles.txt   /path/to/sample_name_amplicon1_graph.txt`

Other arguments
- `--ref [hg19, GRCh37 or GRCh38]`: (Required) Choose reference genome version used in generating AA output.
- `-v/--version`: Print version and exit.
- `-o`: Output filename prefix. Default is prefix of `-i` or `-c`.
- `--add_chr_tag`: If you have a mix of hg19 and GRCh37 amplicons, you can set `--ref hg19` and `--add_chr_tag` to classify them all together.
- `--min_cn_flow`: Minumum cycle CN flow to consider as an amplification (default=1).
- `--min_size`: Minimum cycle size (in bp) to consider as valid amplicon (default=5000).
- `--force`: Disable No amp/Invalid class, if possible. Use only when extremely large CN seeds were used in AA amplicon generation.
- `--plotStyle [noplot, individual]`: Produce a radar-style plot of classification strenghts. Default `noplot`. 
