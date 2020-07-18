# Exploring RNA secondary structure representations

## Introduction
Structural characterization of RNAs is a dynamic field, exposing many modelling possibilities. Every model is usually characterized by an encoding in which to include structural information of a molecule ranging from string representations to graphs. Introducing a re-interpretation of the Shannon Information applied on RNA alignments, we propose a new scoring metric, the **Relative Information Gain**, available for any position in an alignment, showing how different levels of detail encoded in the RNA representation contribute differently to expose structural information.

## How to Build PSSM and Blocks

1. How to make a PSSM from a new MBR

* Download the script from the folder script/
* Choose the MBR version (ex. "Zbear_62")
* Select the MBR file from the folder MBR (ex. MBR/MBR_Zbear_62.tsv)
* Select the file with the alphabet mapping from the folder alphabets (ex. alphabets/Zbear.tsv)
* Be sure that the file gapped_fam_dict.pickle is in the same folders of make_PSSM.py script. You can find it in the script/ folder

Now you can run the script with the command.

`python3 make_PSSM.py $mbrVersion $MBR_matrix.tsv $ALPHAMAP.tsv`




2. How to Build Blocks and MBR from RFAM families (sequence and structure)

* Download the script from the folder script/
* Download the $sequence_folder with sequences as in folder: seq_str_families/sequence/ 
* Download the $structure_folder with sequences in bear format as in folder seq_str_families/bear/
* Download the RFAM seed sequences ( Rfam_no_double.seed ) from the main folder
* Choose the $identity_score threshold (ex. 62)
* Choose the $seq_threshold as the minimum number of RNAs in a RFAM family (ex. 5)
* Choose the $alphabet as in folder alphabets/ (ex. alphabets/Zbear.tsv)
* Indicate the $info_file as output file that will collect all the informations of the built blocks (ex info_Zbear_62.txt)

`python3 make_MBR.py $sequence_folder $structure_folder Rfam_no_double.seed $identity_score $seq_threshold $alphabet $info_file`


