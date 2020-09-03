# Exploring RNA secondary structure representations

## Introduction
Structural characterization of RNAs is a dynamic field, offering many modelling possibilities. Every model is usually 
characterized by an encoding that depicts structural information of the molecule through string representations or graphs. 
Introducing a re-interpretation of the Shannon Information applied on RNA alignments, we propose a new scoring metric, 
the **Relative Information Gain** (RIG). The RIG score is available for any position in an alignment, showing how different 
levels of detail encoded in the RNA representation can contribute differently to convey structural information.

Computed RIG scores presented in this work can be found in the [RIG](data/RIG) folder.

## Installation

### 1. Prepare your system

TODO: add instructions to satisfy all scripts' dependencies

### 2. Clone and enter the repository

```
git clone https://github.com/citterich-lab/RIG.git
cd RIG
```

## Usage


### Plots generation
 
#### - RIG scores with WUSS notation from secondary structure consensus
Execute

```
python3 scripts/plot_RIG_with_WUSS_notation.py 
```

The plots will be generated in the [RIG_WUSS](plots/RIG_WUSS) folder.


#### - RIG scores together with (rescaled) sequence Entropy

Execute

```
python3 scripts/XXXXX.py 
```

to compute the (rescaled) entropy values.

Execute

```
python3 scripts/XXXXX.py 
```

The plots will be generated in the [RIG_Entropy](plots/RIG_Entropy) folder.


#### - RIG scores together with R-scape power values
Execute

```
python3 scripts/plot_RIG_and_RscapePower.py 
```

The plots will be generated in the [RIG_RscapePower](plots/RIG_RscapePower) folder.


### Build RNA Blocks from Rfam families alignment

To build RNA Blocks from structural alignments you need to specify:

* the `structural_alignment_folder` as in folder [bear_alignment](data/bear_alignment) 
* the `identity_score threshold`
* the `alphabet` as in the [alphabets](data/alphabets) folder

#### Example:
```
python3 scripts/make_RNA_Blocks.py data/bear_alignment/bear_new_alignment_Zbear_90/ 90 alphabets/Zbear.tsv
```

The RNA Blocks will be built in the [RNA_Blocks](outputs/RNA_Blocks) folder.


### Build Matrix of Bear encoded RNA (MBR) from Rfam Blocks

To build a MBR you need to specify:

* the `RNA_blocks_folder` (built in the previous **RNA Blocks from Rfam families alignment** step)
* the `identity_score threshold`
* the `alphabet` as in the [alphabets](data/alphabets) folder
* the `info_file` as output file that will collect all the information of the built blocks

#### Example:
```
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_Zbear_90/ 90 data/alphabets/Zbear.tsv Zbear_90
```

The MBR will be built in the [MBRs](outputs/MBRs) folder.


### Build a structural Position Specific Scoring Matrix (sPSSM) from a new Matrix of Bear encoded RNA (MBR)

To build a sPSSM you need to specify:

* the MBR version (for example, `Zbear_90`)
* the MBR file (it can be one of the generated MBRs in the [MBRs](outputs/MBRs) folder)
* the `alphabet` as in the [alphabets](data/alphabets) folder 
* the [gapped_fam_dict.pickle](scripts/gapped_fam_dict.pickle)

#### Example:
```
python3 scripts/make_PSSM.py Zbear_90 outputs/MBRs/MBR_Zbear_90.tsv data/alphabets/Zbear.tsv scripts/gapped_fam_dict.pickle
```

The sPSSM will be built in the [sPSSM](outputs/sPSSMs) folder.




### BlastClust and structural alignment from Rfam families alignment

If you prefer not to use our precomputed structural alignments, you can perform it as follows.

#### Dependencies: blastclust

**Note**: the following instructions are for the **Linux** operating system. Please change the ftp link according to 
your operating system to download the correct version of blastclust 
(see ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/).

```
cd ~
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-ia32-linux.tar.gz
tar -xvf blast-2.2.26-ia32-linux.tar.gz
```

To calculate the alignments, you need to specify:

* the `sequence_folder` as in the [sequence](seq_str_families/sequence) folder
* the `structure_folder` with sequences in BEAR format as in the [bear](seq_str_families/bear) folder
* the Rfam seed sequences (`Rfam_no_double.seed` in the main folder)
* the `identity_score threshold`
* the `seq_threshold` as the minimum number of RNAs in a Rfam family
* the `alphabet` as in the [alphabets](data/alphabets) folder
* the `bin`directory where it is installed `blastclust`

#### Example:
```
cd scripts
python3 scripts/BlustClust_filter_alignment.py seq_str_families/sequence/ seq_str_families/bear/ Rfam_no_double.seed 90 5 data/alphabets/Zbear.tsv ~/blast-2.2.26/bin/
```


### Convert bear files to other alphabets
To convert a file from `fastB` format ([Mattei et al., 2015](https://academic.oup.com/nar/article/43/W1/W493/2467934)) 
to other structural alphabets, execute

```
python3 scripts/mapping.py Zbear.tsv file_bear.fastB > out_file.fastB
```