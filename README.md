# Exploring RNA secondary structure representations

## Introduction
Structural characterization of RNAs is a dynamic field, offering many modelling possibilities. Every model is usually 
characterized by an encoding that depicts structural information of the molecule through string representations or graphs. 
Introducing a re-interpretation of the Shannon Information applied on RNA alignments, we propose a new scoring metric, 
the **Relative Information Gain** (RIG). The RIG score is available for any position in an alignment, showing how different 
levels of detail encoded in the RNA representation can contribute differently to convey structural information.

Computed RIG scores presented in this work can be found in the [RIG](data/Rfam13.0_old/RIG) folder.

## Installation

### 1. Prepare your system

``` pip install -r requirements.txt```

### 2. Clone and enter the repository

```
git clone https://github.com/citterich-lab/RIG.git
cd RIG
```

## Usage

### Short guide

To calculate all RIG scores for all the `alphabet` as in the [alphabets](data/alphabets) folder, removing redundant
primary sequences up to a 50% and 90% of identity, execute 
 ```
bash scripts/calculate_rig_scores.sh
```

#### Step by step guide
##### Build RNA Blocks from Rfam families alignment

To build RNA Blocks from structural alignments you need to specify:

* the `structural_alignment_folder` as in folder [bear_alignment](data/Rfam13.0/bear_alignment) 
* the `identity_score_threshold`
* the `alphabet` as in the [alphabets](data/alphabets) folder

###### Examples:
```
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_bear_90 90 alphabets/bear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_qbear_90 90 alphabets/qbear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_Zbear_90 90 alphabets/Zbear.tsv

python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_bear_62 62 alphabets/bear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_qbear_62 62 alphabets/qbear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_Zbear_62 62 alphabets/Zbear.tsv

python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_bear_50 50 alphabets/bear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_qbear_50 50 alphabets/qbear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_Zbear_50 50 alphabets/Zbear.tsv
```

The RNA Blocks will be built in the [RNA_Blocks](outputs/RNA_Blocks) folder.


##### Build Matrix of Bear encoded RNA (MBR) from Rfam Blocks

To build a MBR you need to specify:

* the `RNA_blocks_folder` (built in the previous **RNA Blocks from Rfam families alignment** step)
* the `identity_score_threshold`
* the `alphabet` as in the [alphabets](data/alphabets) folder
* the `info_file` as output file that will collect all the information of the built blocks

###### Examples:
```
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_bear_90 90 data/alphabets/bear.tsv bear_90
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_qbear_90 90 data/alphabets/qbear.tsv qbear_90
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_Zbear_90 90 data/alphabets/Zbear.tsv Zbear_90

python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_bear_62 62 data/alphabets/bear.tsv bear_62
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_qbear_62 62 data/alphabets/qbear.tsv qbear_62
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_Zbear_62 62 data/alphabets/Zbear.tsv Zbear_62

python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_bear_50 50 data/alphabets/bear.tsv bear_50
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_qbear_50 50 data/alphabets/qbear.tsv qbear_50
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_Zbear_50 50 data/alphabets/Zbear.tsv Zbear_50
```

The MBR will be built in the [MBRs](outputs/MBRs) folder.


##### Build a structural Position Specific Scoring Matrix (sPSSM) from a new Matrix of Bear encoded RNA (MBR)

To build a sPSSM you need to specify:

* the MBR version (for example, `Zbear_90`)
* the MBR file (it can be one of the generated MBRs in the [MBRs](outputs/MBRs) folder)
* the `alphabet` as in the [alphabets](data/alphabets) folder 
* the [gapped_fam_dict.pickle.gz](data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz) file

##### Examples:
```
python3 scripts/make_PSSM.py bear_90 outputs/MBRs/bear_90/MBR_bear_90.tsv data/alphabets/bear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py qbear_90 outputs/MBRs/qbear_90/MBR_qbear_90.tsv data/alphabets/qbear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py Zbear_90 outputs/MBRs/Zbear_90/MBR_Zbear_90.tsv data/alphabets/Zbear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz

python3 scripts/make_PSSM.py bear_62 outputs/MBRs/bear_62/MBR_bear_62.tsv data/alphabets/bear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py qbear_62 outputs/MBRs/qbear_62/MBR_qbear_62.tsv data/alphabets/qbear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py Zbear_62 outputs/MBRs/Zbear_62/MBR_Zbear_62.tsv data/alphabets/Zbear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz

python3 scripts/make_PSSM.py bear_50 outputs/MBRs/bear_50/MBR_bear_50.tsv data/alphabets/bear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py qbear_50 outputs/MBRs/qbear_50/MBR_qbear_50.tsv data/alphabets/qbear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py Zbear_50 outputs/MBRs/Zbear_50/MBR_Zbear_50.tsv data/alphabets/Zbear.tsv data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
```

The sPSSM will be built in the [sPSSM](outputs/sPSSMs) folder.


##### Calculate RIG scores
To calculate the RIG scores, you need to specify:

* the `sPSSM` file (it can be one of the generated MBRs in the [sPSSM](outputs/sPSSM) folder)

##### Examples:
```
python3 scripts/compute_RIG.py outputs/sPSSMs/bear_90/rfam_PSSM_dic_bear_90.pickle
python3 scripts/compute_RIG.py outputs/sPSSMs/qbear_90/rfam_PSSM_dic_qbear_90.pickle
python3 scripts/compute_RIG.py outputs/sPSSMs/Zbear_90/rfam_PSSM_dic_Zbear_90.pickle

python3 scripts/compute_RIG.py outputs/sPSSMs/bear_62/rfam_PSSM_dic_bear_62.pickle
python3 scripts/compute_RIG.py outputs/sPSSMs/qbear_62/rfam_PSSM_dic_qbear_62.pickle
python3 scripts/compute_RIG.py outputs/sPSSMs/Zbear_62/rfam_PSSM_dic_Zbear_62.pickle

python3 scripts/compute_RIG.py outputs/sPSSMs/bear_50/rfam_PSSM_dic_bear_50.pickle
python3 scripts/compute_RIG.py outputs/sPSSMs/qbear_50/rfam_PSSM_dic_qbear_50.pickle
python3 scripts/compute_RIG.py outputs/sPSSMs/Zbear_50/rfam_PSSM_dic_Zbear_50.pickle
```

The RIG scores will be written in the [RIGs](outputs/RIGs) folder.


#### Compute the (normalized) sequence entropy
Execute

```
python3 scripts/compute_entropy.py data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz
```

The output will be written in the [entropy](outputs/entropy) folder.


### Plots generation
 
#### RIG scores with WUSS notation from secondary structure consensus
Execute

```
python3 scripts/plot_RIG_with_WUSS_notation.py
```

The plots will be generated in the [RIG_WUSS](plots/RIG_WUSS) folder.


#### RIG scores minus the (rescaled) sequence entropy
Execute

```
python3 scripts/plot_RIG_minus_Entropy.py
```

The plots will be generated in the [RIG_Entropy](plots/RIG_Entropy) folder.


#### RIG scores together with R-scape power values
Execute

```
python3 scripts/plot_RIG_and_RscapePower.py
```

The plots will be generated in the [RIG_RscapePower](plots/RIG_RscapePower) folder.


### Other

#### BlastClust and structural alignment from Rfam families alignment

If you prefer not to use our precomputed structural alignments, you can perform it as follows.

##### Dependencies: blastclust

**Note**: the following instructions are for the **Linux** operating system. Please change the ftp link according to 
your operating system to download the correct version of blastclust 
(see ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/).

```
cd ~
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-ia32-linux.tar.gz
tar -xvf blast-2.2.26-ia32-linux.tar.gz
```

To calculate the alignments, you need to specify:

* the `sequences_folder` as in the [sequences](data/Rfam13.0/sequences) folder
* the `structure_folder` with sequences in BEAR format as in the [bear](data/Rfam13.0/bear) folder
* the [gapped_fam_dict.pickle.gz](data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz) file
* the `identity_score threshold`
* the `seq_threshold` as the minimum number of RNAs in a Rfam family
* the `alphabet` as in the [alphabets](data/alphabets) folder
* the `bin`directory where it is installed `blastclust`

###### Examples:
```
python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 90 5 data/alphabets/bear.tsv ~/blast-2.2.26/bin/
python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 90 5 data/alphabets/qbear.tsv ~/blast-2.2.26/bin/
python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 90 5 data/alphabets/Zbear.tsv ~/blast-2.2.26/bin/

python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 62 5 data/alphabets/bear.tsv ~/blast-2.2.26/bin/
python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 62 5 data/alphabets/qbear.tsv ~/blast-2.2.26/bin/
python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 62 5 data/alphabets/Zbear.tsv ~/blast-2.2.26/bin/

python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 50 5 data/alphabets/bear.tsv ~/blast-2.2.26/bin/
python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 50 5 data/alphabets/qbear.tsv ~/blast-2.2.26/bin/
python3 scripts/BlustClust_filter_alignment.py data/Rfam13.0/sequences/ data/Rfam13.0/bear/ data/Rfam13.0/gapped_fam/gapped_fam_dict.pickle.gz 50 5 data/alphabets/Zbear.tsv ~/blast-2.2.26/bin/
```

The alignment will be generated in the [alignments](outputs/alignments) folder.


#### Convert bear files to other alphabets
To convert a file from `fastB` format ([Mattei et al., 2015](https://academic.oup.com/nar/article/43/W1/W493/2467934)) 
to other structural alphabets, execute

```
python3 scripts/mapping.py Zbear.tsv file_bear.fastB > out_file.fastB
```