import os

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import seaborn as sns
import sys

import function_mbr

sns.set(context='paper', style='whitegrid', palette='deep', font='serif', font_scale=3, color_codes=True, rc=None)

# Paths and directories
WUSS_path = sys.argv[1]

RIG_dir = "outputs/RIGs/"

bear90_path = os.path.join(RIG_dir, "bear_90_RIGs.tsv")
bear50_path = os.path.join(RIG_dir, "bear_50_RIGs.tsv")
qbear90_path = os.path.join(RIG_dir, "qbear_90_RIGs.tsv")
qbear50_path = os.path.join(RIG_dir, "qbear_50_RIGs.tsv")
zbear90_path = os.path.join(RIG_dir, "zbear_90_RIGs.tsv")
zbear50_path = os.path.join(RIG_dir, "zbear_50_RIGs.tsv")


# Colors WUSS dictionary
WUSS_color_dict = function_mbr.load_and_color_WUSS_dictionary(WUSS_path)

# Load RIG values
RIG_dict = function_mbr.load_rig_values([
    [bear90_path, 'bear90'],
    [bear50_path, 'bear50'],
    [qbear90_path, 'qbear90'],
    [qbear50_path, 'qbear50'],
    [zbear90_path, 'zbear90'],
    [zbear50_path, 'zbear50'],
])

# Remove gaps
with open(WUSS_path) as f:
    for line in f:
        RF, WUSS = line.strip('\n').split('\t')

        for encoding, rf_to_rigs_dict in RIG_dict.items():
            if RF in rf_to_rigs_dict:
                RIG_dict[encoding][RF] = function_mbr.removeGapsFromRIG(WUSS, rf_to_rigs_dict[RF])
            else:
                print('Missing', RF)

##Dati di RIG - entropy RIG
"""
bearRIG-- conservazione generale in struttura, anche solo in senso di evoluzione
zbearRIG-- conservazione stretta della tipologia di elemento
"""

def plot_RF_delta_RIG_ENT(RF, RIG_dict, XXX_dict, WUSS_color_dict, encodings, XXX, filename='test'):
    for encoding in encodings:
        plt.clf()

        if len(XXX_dict[RF]) != len(RIG_dict[encoding][RF]):
            print('None')
            return None

        fig, ax = plt.subplots(1, figsize=(20, 10), sharex='all')
        s = np.array([x - y for x, y in zip(RIG_dict[encoding][RF], XXX_dict[RF])])

        sns.heatmap([s], ax=ax,
                    xticklabels=len(s) // 10 if len(s) > 20 else len(s),
                    yticklabels=False,
                    vmin=-1, vmax=1, cmap='RdBu_r')

        # Generate the subdivision in colors of the various blocks
        colored_chunks = function_mbr.get_colored_chunks(WUSS_color_dict[RF])

        # For each color block, color under the corresponding area
        for count, colored_chunk in enumerate(colored_chunks):
            plt.fill_between(colored_chunk[0], 0.75, y2=1, facecolor=colored_chunk[2], alpha=1)
        plt.axhline(0.75, color='k')

        ax.set_title('{}'.format(RF), fontsize=28)
        ax.set_ylabel('{} RIG - {} RIG'.format(encoding, XXX), fontsize=28)
        ax.set_xlabel('position in RF alignment', fontsize=28)

        plt.xticks(rotation=45)
        # ax.set_ylim([-1, 1])

        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)

        # Add the legend
        legend_elements = [Patch(facecolor='r', label='Stem'),
                           Patch(facecolor='g', label='Hairpin Loop'),
                           Patch(facecolor='c', label='Bulge/Interior'),
                           Patch(facecolor='b', label='Others'),
                           Patch(facecolor='m', label='Pseudoknot')]

        legend = ax.legend(handles=legend_elements, frameon=1, fontsize=22)
        frame = legend.get_frame()
        frame.set_facecolor('w')

        plt.tight_layout()
        # plt.savefig(filename + f"_{encoding}.pdf")
        plt.savefig(filename + f"_{encoding}_{XXX}.png", ppi=600)

        plt.close()


if sys.argv[2] == 'entropy':
    path_entropy = 'outputs/entropy/entropy.tsv'

    ENT_dict = {}
    with open(path_entropy) as f:
        for line in f:
            # RF\tENT
            dat = line.split()
            ENT_dict[dat[0]] = [x for x in map(float, dat[1:])]


    RIG_vs_ENT_output_dir = 'outputs/plots/RIG_Entropy/'
    if not os.path.exists(RIG_vs_ENT_output_dir):
        os.makedirs(RIG_vs_ENT_output_dir)


    for RF_ in RIG_dict[list(RIG_dict.keys())[0]]:
        print(RF_)

        plot_RF_delta_RIG_ENT(
            RF_,
            RIG_dict=RIG_dict,
            XXX_dict=ENT_dict,
            WUSS_color_dict=WUSS_color_dict,
            encodings=['bear90', 'qbear90', 'zbear90'],
            XXX='entropy',
            filename=f"{RIG_vs_ENT_output_dir}{RF_}_90"
        )

else:
    RIG_vs_RIG_output_dir = 'outputs/plots/RIG_RIG/'
    if not os.path.exists(RIG_vs_RIG_output_dir):
        os.makedirs(RIG_vs_RIG_output_dir)

    for RF_ in RIG_dict[list(RIG_dict.keys())[0]]:
        print(RF_)

        plot_RF_delta_RIG_ENT(
            RF_,
            RIG_dict=RIG_dict,
            XXX_dict=RIG_dict['zbear90'],
            WUSS_color_dict=WUSS_color_dict,
            encodings=['bear90'],
            XXX='zbear90',
            filename=f"{RIG_vs_RIG_output_dir}{RF_}_90"
        )
