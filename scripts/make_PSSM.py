import os
import sys
import gzip
import argparse
import pickle
import pandas as pd

from matrix2 import buildPSSM_alphabet

dir_output_sPSSMs = 'outputs/sPSSMs'

def mapAligns(gapped_fam_dict, alphamap):
    """convert alignments in the new alphabet"""
    for RF in gapped_fam_dict:
        # print(RF)
        for ID in gapped_fam_dict[RF]:
            if gapped_fam_dict[RF][ID].get('bear'):
                gapped_fam_dict[RF][ID]['alpha'] = "".join(
                    [alphamap[ch] for ch in gapped_fam_dict.get(RF).get(ID).get('bear')]
                )


mbrVersion = sys.argv[1]

MBR = sys.argv[2]

ALPHAMAP = sys.argv[3]

GAPFAMDICT = sys.argv[4]

ignore_these_families = [
    #'RF00210', 'RF01879', #RF is missing
    #'RF02767', 'RF02768', 'RF02770', 'RF02773', 'RF02775', 'RF02781', 'RF02783'
]

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mbr", help="the substitution matrix to test")
parser.add_argument("-a", "--alpha", help="pickle of alphabet dictionary (must correspond to the alphabet \
                                            used with the substitution matrix). With respect to standard BEAR")
parser.add_argument("-gfd", "--gapfamdict", help="the rfam alignments with gaps and bear pickle (gzipped)")
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

args = parser.parse_args(args=['--alpha', ALPHAMAP, '-v', '-m', MBR, '-gfd', GAPFAMDICT])

with gzip.open(args.gapfamdict, 'rb') as afile:
    gapped_fam_dict = pickle.load(afile)

# Prepare the alphabets, adding the gap symbol
f = open(args.alpha).readlines()
alph_bear = {x.split()[0]: x.split()[1] for x in f}
alph_bear['-'] = '-'


alphamaps = [{key: alph_bear[KEY] for key in KEY} for KEY in alph_bear]
alphamap = {}
for ap in alphamaps:
    alphamap.update(ap)


mbr = pd.read_table(args.mbr, index_col=0)

mapAligns(gapped_fam_dict, alphamap)

rfams = {}

num_families_div_10 = len(gapped_fam_dict) // 10

for num_fam, RF in enumerate(gapped_fam_dict):
    # if num_fam % num_families_div_10 == 0:
    #    print('{}% ({}/{})'.format((num_fam // num_families_div_10) * 10, num_fam, len(gapped_fam_dict)))

    if RF not in ignore_these_families:
        rfams[RF] = [
            [gapped_fam_dict[RF][ID]['sequence'] for ID in gapped_fam_dict[RF]],
            [gapped_fam_dict[RF][ID]['alpha'] for ID in gapped_fam_dict[RF]]
        ]


PSSMs_alpha = []
rfam_list = []

num_families_div_20 = len(rfams) // 20

# use RFAMS
print('Build PSSM: STARTED')
for num_fam, rfam in enumerate(rfams):
    if num_fam % num_families_div_20 == 0:
        print('{:.2f}% ({}/{})'.format((num_fam / len(rfams)) * 100.0, num_fam, len(rfams)))

    rfam_list.append(rfam)
    PSSM_b = buildPSSM_alphabet(rfams[rfam], mbr)

    PSSMs_alpha.append(PSSM_b)
print('Build PSSM: DONE!')


dir_output_sPSSMs_name_id = os.path.join(dir_output_sPSSMs, mbrVersion)

if not os.path.exists(dir_output_sPSSMs_name_id):
    os.makedirs(dir_output_sPSSMs_name_id)

with gzip.open(os.path.join(dir_output_sPSSMs_name_id, 'rfam_PSSM_{}.pickle.gz'.format(mbrVersion)), 'wb') as handle:
    pickle.dump(PSSMs_alpha, handle)

dic_PSSM = {}
for i in range(0, len(rfam_list)):
    dic_PSSM[rfam_list[i]] = PSSMs_alpha[i]

with gzip.open(os.path.join(dir_output_sPSSMs_name_id, 'rfam_PSSM_dic_{}.pickle.gz'.format(mbrVersion)), 'wb') as handle:
    pickle.dump(dic_PSSM, handle)
