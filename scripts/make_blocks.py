import function_mbr as f
import sys

#How to run python make_blocks.py bear_new_alignment_Zbear_62/ 62 Zbear.tsv

folder_alignment=sys.argv[1]
#'bear_new_alignment_'+name+'_'+id_blustClust+'/'
id_blustClust=sys.argv[2]
file_alph=sys.argv[3]

name=file_alph.split('.')[0]
f.make_blocks(folder_alignment, name, id_blustClust)