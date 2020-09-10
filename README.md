## AmpliconClassifier

#### Classify AmpliconArchitect amplicons to determine what type of focal amplification is present.

**Current version: 0.2.5**

**1. Prerequisites**: 
- Supports both python2 and python3.
- intervaltree (python library):  `pip install intervaltree`
- `$AA_DATA_REPO` environment variable. For instructions see [AmpliconArchitect installation](https://github.com/jluebeck/AmpliconArchitect#data-repositories). 
- matplotlib (python library, optional): `pip install matplotlib`

**2. Usage**:

To generate classifications for a list of amplicons, I recommend doing

`python amplicon_classifier.py --ref [hg19, GRCh37, or GRCh38] --input [file with list of your amplicons] > classifier_stdout.log`

the `--input` argument must be formatted as follows, with one amplicon per line:

`sample_name_amplicon1   /path/to/sample_name_amplicon1_cycles.txt   /path/to/sample_name_amplicon1_graph.txt`

Alternatively, you can just classify one amplicon at a time

`python amplicon_classifier.py --ref [hg19, GRCh37, or GRCh38] --cycles [/path/to/amplicon_cycles.txt] --graph [/path/to/amplicon_graph.txt] > classifier_stdout.log`

There is also an experimental option you can set to visualize the strength of each amplicon class assigned to an amplicon, which can be turned on by setting `--plotStyle individual`.

If you have data from both GRCh37 and hg19 that you are combining, you can set the flag `--add_chr_tag` to add the "chr" prefix to each chromosome name and effectively unify everything as hg19-based.

**3. Output**:

The most important file will be the file `[output_prefix]_amplicon_classification_profiles.tsv`. 
This contains an abstract classification of the amplicon, and also indicates in separate columns "BFB+" and "ecDNA+" status.

**4. Command line arguments**:

Running on a single AA amplicon:
- `-c/--cycles`: AA cycles file
- `-g/--graph`: AA graph file

OR running on multiple amplicons
- `-i/--input`: Tab-separated file containing one or more amplicons formatted as

`sample_name_amplicon1   /path/to/sample_name_amplicon1_cycles.txt   /path/to/sample_name_amplicon1_graph.txt`

Other arguments
- `--ref [hg19, GRCh37 or GRCh38]`: (Required) Choose reference genome version used in generating AA output.
- `-v/--version`: Print version ane exit.
- `-o`: Output filename prefix. Default is prefix of `-i` or `-c`.
- `--add_chr_tag`: If you have a mix of hg19 and GRCh37 amplicons, you can set `--ref hg19` and `--add_chr_tag` to classify them all together.
- `--min_cn_flow`: Minumum cycle CN flow to consider as an amplification (default=1).
- `--min_size`: Minimum cycle size (in bp) to consider as valid amplicon (default=5000).
- `--force`: Disable No amp/Invalid class, if possible. Use only when extremely large CN seeds were used in AA amplicon generation.
- `--plotStyle [noplot, individual]`: Produce a radar-style plot of classification strenghts. Default `noplot`. 