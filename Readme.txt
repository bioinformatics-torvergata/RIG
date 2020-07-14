############################### FOLDERS & FILES ###############################

#folders with bear alignment based on alphabets and identity score
#bear_alignment
	|--bear_alignment___62
	|--bear_alignment___50
	|--bear_alignment___90
	|
	|--bear_new_alignment_bear_50
	|--bear_new_alignment_bear_62
	|--bear_new_alignment_bear_90
	|
	|--bear_new_alignment_Zbear_50
	|--bear_new_alignment_Zbear_62
	|--bear_new_alignment_Zbear_90
	|
	|--bear_new_alignment_qbear_50
	|--bear_new_alignment_qbear_62
	|--bear_new_alignment_qbear_90


#folders with Blocks based on alphabets and identity score
#----Blocks
	|--blocks_new_bear_bear_50 
	|--blocks_new_bear_bear_62  
	|--blocks_new_bear_bear_90       
	|
	|--blocks_new_bear_Zbear_50 
	|--blocks_new_bear_Zbear_62  
	|--blocks_new_bear_Zbear_90   
	|
	|--blocks_new_bear_qbear_50 
	|--blocks_new_bear_qbear_62  
	|--blocks_new_bear_qbear_90           



#folder with all the RNA sequences/structure(bear) --> one file for each RFAM family
#seq_str_families
	|--sequence
	|--bear (structure)


######## Rfam_no_double.seed ##########
#seed file downloaded from Rfam


######## Table with new characters conversion
#alphabets
       	|----bear.tsv
	|----Zbear.tsv
	|----Qbear.tsv
	|----readme (how to build a customized structural alphabet)



######## MBR from all the alphabets and different score identity
## MBR
   |---MBR.tsv
   |---MBR_Zbear_50.tsv
   |---MBR_Zbear_62.tsv
   |---MBR_Zbear_90.tsv
   |---MBR_qbear_50.tsv
   |---MBR_qbear_62.tsv
   |---MBR_qbear_90.tsv
   |---MBR_bear_50.tsv
   |---MBR_bear_62.tsv
   |---MBR_bear_90.tsv



#PSSM built for all the alphabets in pickle file format 
#--PSSM as dict
	|--rfam_PSSM_dic_bear_50.pickle
	|--rfam_PSSM_dic_bear_62.pickle
	|--rfam_PSSM_dic_bear_90.pickle
	|
	|--rfam_PSSM_dic_Zbear_50.pickle
	|--rfam_PSSM_dic_Zbear_62.pickle
	|--rfam_PSSM_dic_Zbear_90.pickle
	|
	|--rfam_PSSM_dic_qbear_50.pickle
	|--rfam_PSSM_dic_qbear_62.pickle
	|--rfam_PSSM_dic_qbear_90.pickle

#--PSSM as matrix
	|--rfam_PSSM_bear_50.pickle
	|--rfam_PSSM_bear_62.pickle
	|--rfam_PSSM_bear_90.pickle
	|
	|--rfam_PSSM_Zbear_50.pickle
	|--rfam_PSSM_Zbear_62.pickle
	|--rfam_PSSM_Zbear_90.pickle
	|
	|--rfam_PSSM_qbear_50.pickle
	|--rfam_PSSM_qbear_62.pickle
	|--rfam_PSSM_qbear_90.pickle



################# SCRIPT ##########################

#Make PSSM from a new MBR (MARCO)
python3 make_PSSM.py $mbrVersion $MBR_matrix.tsv $ALPHAMAP.tsv (need gapped_fam_dict.pickle in same folder)

#Build Blocks and MBR from RFAM families (sequence and structure)
#Make_MBR.py
#How to run (Example): python make_MBR.py seq_str_families/sequence/ seq_str_families/bear/ Rfam_no_double.seed 90 5 alphabets/bear.tsv info_bear_90.txt

