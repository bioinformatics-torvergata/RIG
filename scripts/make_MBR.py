import sys

import function_mbr as f

blocks_folder = sys.argv[1]
id_blustClust = sys.argv[2]
file_alph = sys.argv[3]
file_info = sys.argv[4]

f.Make_MBR_from_blocks(blocks_folder, id_blustClust, file_alph, file_info)
