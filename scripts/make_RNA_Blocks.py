import os
import sys

import function_mbr as f

folder_alignment = sys.argv[1]
id_blustClust = sys.argv[2]
file_alph = sys.argv[3]

name = os.path.basename(file_alph).split('.')[0]

f.make_blocks(folder_alignment, name, id_blustClust)
