# RIG

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
* Download the sequences files as in folder: seq_str_families/sequence/
* Download the structures files in bear format as in folder seq_str_families/bear/
* Download the RFAM seed sequences ( Rfam_no_double.seed ) from the main folder
* Choose the identity score threshold (ex. 62)
* Choose the minimum number of RNAs in a RFAM family (ex. 5)
* Choose the alphabet as in folder alphabets/ (ex. alphabets/Zbear.tsv)
* Indicate the file name as output file that will collect all the informations of the built blocks (ex info_Zbear_62.txt)

`python make_MBR.py seq_str_families/sequence/ seq_str_families/bear/ Rfam_no_double.seed 62 5 alphabets/Zbear.tsv info_Zbear_62.txt`


