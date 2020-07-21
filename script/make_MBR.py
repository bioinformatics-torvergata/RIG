import function_mbr as f
import sys

#How_to_run
#python make_MBR.py blocks_new_bear_Zbear_62/ 62 Zbear.tsv prova.txt
    
blocks_folder=sys.argv[1]
id_blustClust=sys.argv[2]
file_alph=sys.argv[3]
file_info=sys.argv[4]


f.Make_MBR_from_blocks(blocks_folder, id_blustClust, file_alph, file_info)
