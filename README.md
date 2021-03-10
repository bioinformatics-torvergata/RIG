# Relative Information Gain: Shannon entropy-based measure of the relative structural conservation in RNA alignments

## Introduction
Structural characterization of RNAs is a dynamic field, offering many modelling possibilities. RNA secondary structure
models are usually characterized by an encoding that depicts structural information of the molecule through string
representations or graphs. Introducing a re-interpretation of the Shannon Information applied on RNA alignments, we propose
a new scoring metric, the **Relative Information Gain** (RIG). The RIG score is available for any position in an alignment,
showing how different levels of detail encoded in the RNA representation can contribute differently to convey structural
information.

Computed RIG scores can be found in the [RIG](outputs/RIGs) folder.

## Citation

**Marco Pietrosanto\*, Marta Adinolfi\*, Andrea Guarracino\*, Fabrizio Ferrè, Gabriele Ausiello, Ilio Vitale, Manuela Helmer-Citterich**. [Relative Information Gain: Shannon Entropy-Based Measure of the Relative Structural Conservation in RNA Alignments](https://doi.org/10.1093/nargab/lqab007), NAR Genomics and Bioinformatics, 2021
**\*Shared first authorship**

## Installation

### 1. Clone and enter the repository

```
git clone https://github.com/citterich-lab/RIG.git
cd RIG
```

### 2. Prepare your system

``` pip3 install -r requirements.txt```

## Usage

### Short guide

With the following instruction the RIG scores are calculated by applying the zBEAR alphabet and removing redundant
primary sequences up to 90% of identity. The computation are executed using a precalculated structural Position Specific 
Scoring Matrix (sPSSM) [available](outputs/sPSSMs/zbear_90/rfam_PSSM_dic_zbear_90.pickle.gz) in the repository.

**Note**: for a complete guide on how to calculate the RIG scores from scratch by applying all the encodings in the
[alphabets](data/alphabets) folder and removing redundant primary sequences up to 50%, 62%, and 90% of identity, refer
to the next paragraph [Complete guide](#complete-guide).

To calculate the RIG scores, you just need to specify the `sPSSM` file (it can be one of the generated matrices in the
[sPSSM](outputs/sPSSM) folder).

```
python3 scripts/compute_RIG.py outputs/sPSSMs/zbear_90/rfam_PSSM_dic_zbear_90.pickle.gz
```

The RIG scores will be written in the [zbear_90_RIGs.tsv](outputs/RIGs/zbear_90_RIGs.tsv) file.


### [Complete guide](#complete-guide])

In the following there are instructions to calculate the RIG scores from scratch by applying all the encodings in the
[alphabets](data/alphabets) folder and removing redundant primary sequences up to 50%, 62%, and 90% of identity.


#### Preparation to a new Rfam release

If needed, to create the data for a new Rfam release you need to specify:

* the new gzipped Rfam seed (the last one was [Rfam.seed.gz](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.2/Rfam.seed.gz));
* the output directory;
* the path of the RNAfold program (you can find an executable for Linux operative systems in the [tools](tools) folder);
* the path of the BearEncoder program (you can find the jar in the [tools](tools) folder).

```
python3 scripts/prepare_new_rfam_release.py directory/new/gzipped/rfam-seed/Rfam.seed.gz data/Rfam14.2/ tools/RNAfold tools/BearEncoder.jar
```

#### Structural alignment from Rfam families alignment

##### Dependencies: blastclust

**Note**: the following instructions are for the **Linux** operating system. Please change the ftp link according to 
your operating system to download the correct version of blastclust 
(see ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/).

```
cd ~
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-ia32-linux.tar.gz
tar -xvf blast-2.2.26-ia32-linux.tar.gz
```

To create the BEAR alignments, you need to specify:

* the [sequences](data/Rfam14.2/sequences) folder;
* the [bear_alignment](data/Rfam14.2/bear_alignment) folder with the structural alignments encoded in BEAR;
* the [gapped_fam_dict.pickle.gz](data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz) file;
* the identity threshold;
* the `seq_threshold` as the minimum number of RNAs in a Rfam family;
* the `alphabet` as in the [alphabets](data/alphabets) folder;
* the path of the `blastclust` program.

```
python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 90 5 data/alphabets/bear.tsv ~/blast-2.2.26/bin/blastclust
python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 90 5 data/alphabets/qbear.tsv ~/blast-2.2.26/bin/blastclust
python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 90 5 data/alphabets/zbear.tsv ~/blast-2.2.26/bin/blastclust

python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 62 5 data/alphabets/bear.tsv ~/blast-2.2.26/bin/blastclust
python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 62 5 data/alphabets/qbear.tsv ~/blast-2.2.26/bin/blastclust
python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 62 5 data/alphabets/zbear.tsv ~/blast-2.2.26/bin/blastclust

python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 50 5 data/alphabets/bear.tsv ~/blast-2.2.26/bin/blastclust
python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 50 5 data/alphabets/qbear.tsv ~/blast-2.2.26/bin/blastclust
python3 scripts/BlustClust_filter_alignment.py data/Rfam14.2/sequences/ data/Rfam14.2/bear/ data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz 50 5 data/alphabets/zbear.tsv ~/blast-2.2.26/bin/blastclust
```

The alignment will be generated in the [alignments](outputs/alignments) folder.


#### Build RNA Blocks from Rfam families alignment

To build RNA Blocks from structural alignments you need to specify:

* the [bear_alignment](data/Rfam14.2/bear_alignment) folder with the structural alignments encoded in BEAR;
* the identity threshold;
* the `alphabet` as in the [alphabets](data/alphabets) folder.

```
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_bear_90 90 alphabets/bear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_qbear_90 90 alphabets/qbear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_zbear_90 90 alphabets/zbear.tsv

python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_bear_62 62 alphabets/bear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_qbear_62 62 alphabets/qbear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_zbear_62 62 alphabets/zbear.tsv

python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_bear_50 50 alphabets/bear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_qbear_50 50 alphabets/qbear.tsv
python3 scripts/make_RNA_Blocks.py outputs/alignments/bear_new_alignment_zbear_50 50 alphabets/zbear.tsv
```

The RNA Blocks will be built in the [RNA_Blocks](outputs/RNA_Blocks) folder.


#### Build Matrix of Bear encoded RNA (MBR) from Rfam Blocks

To build a MBR you need to specify:

* the [RNA_Blocks](outputs/RNA_Blocks) folder (built in the previous **RNA Blocks from Rfam families alignment** step)
* the identity threshold;
* the `alphabet` as in the [alphabets](data/alphabets) folder;
* the name of the `info_file` that will collect all the information of the built blocks.

```
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_bear_90 90 data/alphabets/bear.tsv bear_90
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_qbear_90 90 data/alphabets/qbear.tsv qbear_90
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_zbear_90 90 data/alphabets/zbear.tsv zbear_90

python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_bear_62 62 data/alphabets/bear.tsv bear_62
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_qbear_62 62 data/alphabets/qbear.tsv qbear_62
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_zbear_62 62 data/alphabets/zbear.tsv zbear_62

python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_bear_50 50 data/alphabets/bear.tsv bear_50
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_qbear_50 50 data/alphabets/qbear.tsv qbear_50
python3 scripts/make_MBR.py outputs/RNA_Blocks/blocks_new_bear_zbear_50 50 data/alphabets/zbear.tsv zbear_50
```

The MBR will be built in the [MBRs](outputs/MBRs) folder.


#### Build a structural Position Specific Scoring Matrix (sPSSM) from a new Matrix of Bear encoded RNA (MBR)

To build a sPSSM you need to specify:

* the MBR version (for example, `zbear_90`);
* the MBR file (it can be one of the generated MBRs in the [MBRs](outputs/MBRs) folder);
* the `alphabet` as in the [alphabets](data/alphabets) folder;
* the [gapped_fam_dict.pickle.gz](data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz) file.

```
python3 scripts/make_PSSM.py bear_90 outputs/MBRs/bear_90/MBR_bear_90.tsv data/alphabets/bear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py qbear_90 outputs/MBRs/qbear_90/MBR_qbear_90.tsv data/alphabets/qbear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py zbear_90 outputs/MBRs/zbear_90/MBR_zbear_90.tsv data/alphabets/zbear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz

python3 scripts/make_PSSM.py bear_62 outputs/MBRs/bear_62/MBR_bear_62.tsv data/alphabets/bear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py qbear_62 outputs/MBRs/qbear_62/MBR_qbear_62.tsv data/alphabets/qbear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py zbear_62 outputs/MBRs/zbear_62/MBR_zbear_62.tsv data/alphabets/zbear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz

python3 scripts/make_PSSM.py bear_50 outputs/MBRs/bear_50/MBR_bear_50.tsv data/alphabets/bear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py qbear_50 outputs/MBRs/qbear_50/MBR_qbear_50.tsv data/alphabets/qbear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
python3 scripts/make_PSSM.py zbear_50 outputs/MBRs/zbear_50/MBR_zbear_50.tsv data/alphabets/zbear.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
```

The sPSSM will be built in the [sPSSM](outputs/sPSSMs) folder.


#### Calculate RIG scores
To calculate the RIG scores, you need to specify:

* the `sPSSM` file (it can be one of the generated matrices in the [sPSSM](outputs/sPSSM) folder).

```
python3 scripts/compute_RIG.py outputs/sPSSMs/bear_90/rfam_PSSM_dic_bear_90.pickle.gz
python3 scripts/compute_RIG.py outputs/sPSSMs/qbear_90/rfam_PSSM_dic_qbear_90.pickle.gz
python3 scripts/compute_RIG.py outputs/sPSSMs/zbear_90/rfam_PSSM_dic_zbear_90.pickle.gz

python3 scripts/compute_RIG.py outputs/sPSSMs/bear_62/rfam_PSSM_dic_bear_62.pickle.gz
python3 scripts/compute_RIG.py outputs/sPSSMs/qbear_62/rfam_PSSM_dic_qbear_62.pickle.gz
python3 scripts/compute_RIG.py outputs/sPSSMs/zbear_62/rfam_PSSM_dic_zbear_62.pickle.gz

python3 scripts/compute_RIG.py outputs/sPSSMs/bear_50/rfam_PSSM_dic_bear_50.pickle.gz
python3 scripts/compute_RIG.py outputs/sPSSMs/qbear_50/rfam_PSSM_dic_qbear_50.pickle.gz
python3 scripts/compute_RIG.py outputs/sPSSMs/zbear_50/rfam_PSSM_dic_zbear_50.pickle.gz
```

The RIG scores will be written in the [RIGs](outputs/RIGs) folder.


#### Compute the (normalized) sequence entropy
Execute

```
python3 scripts/compute_entropy.py data/Rfam14.2/SS_cons/SS_cons_WUSS.tsv data/Rfam14.2/gapped_fam/gapped_fam_dict.pickle.gz
```

The output will be written in the [entropy](outputs/entropy) folder.


### Plots generation
 
#### RIG scores with WUSS notation from secondary structure consensus
Execute

```
python3 scripts/plot_RIG_with_WUSS_notation.py data/Rfam14.2/SS_cons/SS_cons_WUSS.tsv
```

The plots will be generated in the [RIG_WUSS](outputs/plots/RIG_WUSS) folder.


#### RIG scores minus the (rescaled) sequence entropy
Execute

```
python3 scripts/plot_RIG_minus_EntropyOrRIG.py data/Rfam14.2/SS_cons/SS_cons_WUSS.tsv entropy
```

The plots will be generated in the [RIG_Entropy](outputs/plots/RIG_Entropy) folder.


#### RIG scores minus RIG scores
Execute

```
python3 scripts/plot_RIG_minus_EntropyOrRIG.py data/Rfam14.2/SS_cons/SS_cons_WUSS.tsv RIG
```

The plots will be generated in the [RIG_Entropy](outputs/plots/RIG_Entropy) folder.


#### RIG scores together with R-scape power values

##### Calculate R-scape power

###### Dependencies: R-scape

**Note**: download [here](http://eddylab.org/software/rscape/) the source code distribution of R-scape, and follow the
installation instructions.

```
python3 scripts/compute_rscape_power.py data/Rfam14.2/stockholm/ data/Rfam14.2/Rscape ~/path/where/rscape/executable/is/R-scape
```

Execute

```
python3 scripts/plot_RIG_and_RscapePower.py data/Rfam14.2/SS_cons/SS_cons_WUSS.tsv data/Rfam14.2/Rscape
```

The plots will be generated in the [RIG_RscapePower](outputs/plots/RIG_RscapePower) folder.


### Other

##### Convert bear files to other alphabets
To convert a bear file from `fastB` format ([Mattei et al., 2015](https://academic.oup.com/nar/article/43/W1/W493/2467934)) 
to other structural alphabets, execute

```
python3 scripts/mapping.py data/alphabets/zbear.tsv data/Rfam14.2/bear/RF00001.folded.fastb > RF00001.folded.zbear.fastb
```
