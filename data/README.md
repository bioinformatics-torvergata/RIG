## Data

This folder contains:

- the `MBR.tsv` file, which contains the Matrix of Bear encoded RNA (MBR) from [Mattei et al., 2014](https://academic.oup.com/nar/article/42/10/6146/2436561);
- the [alphabets](alphabets) folder, which contains the tables for mapping all the encodings with the structural elements;
- the [bear_alignment](bear_alignment) folder, which contains the structural alignments expressed with the different 
encodings, and obtained by removing redundant primary sequences up to 50% and up to 90% of similarity for each Rfam seed 
alignment;
- the [Rfam_stockholm_rscapes](Rfam_stockholm_rscapes) folder, which contains the [R-scape](http://eddylab.org/R-scape/) 
outputs. The files were obtained specifying as input the RNA multiple sequence alignments in 
[Stockholm format](https://en.wikipedia.org/wiki/Stockholm_format).
- the [RIG](RIG) folder, which contains the RIG scores calculated using RNA Blocks obtained by removing redundant 
primary sequences up to 50% and up to 90% of similarity for each Rfam seed alignment. In particular,
    - the [RIG/nogaps](RIG/nogaps) folder contains the RIG scores, without the gap positions in the 
    alignments. For each file, the rows are the Rfam seed alignments, the columns are the alignment positions.
    - the [RIG/nogaps](RIG/withgaps) folder contains the RIG scores, including the gap positions in the 
    alignments. For each file, the rows are the Rfam seed alignments, the columns are the alignment positions. Moreover:
        - the [All_RIGs.filled_columns.xlsx](RIG/withgaps/All_RIGs.filled_columns.xlsx) file contains all the RIG
        scores, including the gap positions in the alignment, for both the similarity thresholds, but filling the positions
        with NaN values for shorter alignments. This file exists to simplify the downstream analyses.
- the [SS_cons](SS_cons) folder, which contains the representations of Rfam seed alignment:
    - in dot-bracket notation (file [SS_cons_DB.tsv](SS_cons/SS_cons_DB.tsv)).
    - in WUSS notation (file [SS_cons_WUSS.tsv](SS_cons/SS_cons_WUSS.tsv)).