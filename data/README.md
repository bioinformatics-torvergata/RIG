## Data

This folder contains the files used as input in the scripts.

The `MBR.tsv` file, which contains the Matrix of Bear encoded RNA (MBR) from [Mattei et al., 2014](https://academic.oup.com/nar/article/42/10/6146/2436561);
- the [alphabets](alphabets) folder, which contains the tables for mapping all the encodings with the structural elements;
- the [bear_alignment](Rfam13.0/bear_alignment) folder, which contains the structural alignments;
- the [Rfam_stockholm_rscapes](Rfam13.0/Rfam_stockholm_rscapes) folder, which contains the [R-scape](http://eddylab.org/R-scape/) 
outputs;
- the [RIG](Rfam13.0_old/RIG) folder, which contains the Relative Information Gain (RIG) scores calculated using [RNA Blocks](../outputs/RNA_Blocks)
 obtained by removing redundant primary sequences up to 50% and up to 90% of similarity for each Rfam seed alignment.
- the [SS_cons](Rfam13.0/SS_cons) folder, which contains the representations of Rfam seed alignment in WUSS notation (file [SS_cons_WUSS.tsv](Rfam13.0/SS_cons/SS_cons_WUSS.tsv)).