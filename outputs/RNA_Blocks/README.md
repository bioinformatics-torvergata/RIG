## RNA Blocks

This folder will contain the RNA blocks generated by the [make_RNA_Blocks.py](../../scripts/make_RNA_Blocks.py) script.

For each Rfam seed alignment, we removed redundant primary sequences up to a specified threshold of identity (50%, 62%, 
and 90%), and considered the underlying alignments of secondary structures. Then, we converted the secondary structure 
in each of the structural encoding presented.

The [make_RNA_Blocks.py](../../scripts/make_RNA_Blocks.py) script selects each column in the alignments to be part of the RNA Block 
of that family if
1) no gaps were present, and 
2) a structural consensus, dependent on the chosen alphabet, exists (i.e. there must be a character with a relative
frequency >50%).