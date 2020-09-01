import function_mbr as f
import sys


#How to run
## python BlustClust_filter_alignment.py ../seq_str_families/sequence/ ../seq_str_families/bear/ ../Rfam_no_double.seed 62 5 ../alphabets/Zbear.tsv
folder = sys.argv[1]
folder_bear = sys.argv[2]
RFAM_seed_file = sys.argv[3]
id_blustClust = sys.argv[4]
filter_nSeq = sys.argv[5]
file_alph = sys.argv[6]

f.BlustClust_filter_alignment(folder, folder_bear, RFAM_seed_file, id_blustClust, filter_nSeq, file_alph)