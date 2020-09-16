import os
import sys

dir_input_stockholm = sys.argv[1]
dir_output_rscape = sys.argv[2]
path_rscape = sys.argv[3]

if not os.path.exists(dir_output_rscape):
    os.makedirs(dir_output_rscape)

for xxx_txt in sorted([x for x in os.listdir(dir_input_stockholm) if x .endswith('.txt')]):
    path_xxx_txt = os.path.join(dir_input_stockholm, xxx_txt)
    print(path_xxx_txt)
    os.system(f'{path_rscape} --nofigures -s ' + path_xxx_txt + ' > ' + os.path.join(dir_output_rscape, xxx_txt.replace('.txt', '.rscape.txt')))
