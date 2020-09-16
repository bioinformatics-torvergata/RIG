import sys

import function_mbr as f

folder_seq = sys.argv[1]
folder_bear = sys.argv[2]
path_gapfam_pickle_gz = sys.argv[3]
id_blustClust = sys.argv[4]
filter_nSeq = sys.argv[5]
file_alph = sys.argv[6]
path_blustclust = sys.argv[7]

f.BlustClust_filter_alignment(path_blustclust, folder_seq, folder_bear, path_gapfam_pickle_gz, id_blustClust, filter_nSeq, file_alph)