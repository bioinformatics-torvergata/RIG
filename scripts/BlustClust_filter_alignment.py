import sys

import function_mbr as f

folder_seq = sys.argv[1]
folder_bear = sys.argv[2]
RFAM_seed_file = sys.argv[3]
id_blustClust = sys.argv[4]
filter_nSeq = sys.argv[5]
file_alph = sys.argv[6]
basedir_blustclust = sys.argv[7]

f.BlustClust_filter_alignment(basedir_blustclust, folder_seq, folder_bear, RFAM_seed_file, id_blustClust, filter_nSeq, file_alph)