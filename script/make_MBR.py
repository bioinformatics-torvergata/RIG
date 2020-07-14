import function_mbr as f
import sys

folder = sys.argv[1]
folder_bear = sys.argv[2]
RFAM_seed_file = sys.argv[3]
id_blustClust = sys.argv[4]
filter_nSeq = sys.argv[5]
file_alph = sys.argv[6]
file_info = sys.argv[7]

f.Make_MBR_from_alignment(folder, folder_bear, RFAM_seed_file, id_blustClust, filter_nSeq, file_alph, file_info)
