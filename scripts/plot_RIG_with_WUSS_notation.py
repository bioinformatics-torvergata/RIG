import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import sys

import seaborn as sns
sns.set(context='paper', style='white', palette='deep', font='serif', font_scale=2, color_codes=True, rc=None)

import function_mbr

# Paths and directories
WUSS_path = sys.argv[1]

RIG_dir = "outputs/RIGs/"

bear90_path = os.path.join(RIG_dir, "bear_90_RIGs.tsv")
bear50_path = os.path.join(RIG_dir, "bear_50_RIGs.tsv")
qbear90_path = os.path.join(RIG_dir, "qbear_90_RIGs.tsv")
qbear50_path = os.path.join(RIG_dir, "qbear_50_RIGs.tsv")
zbear90_path = os.path.join(RIG_dir, "zbear_90_RIGs.tsv")
zbear50_path = os.path.join(RIG_dir, "zbear_90_RIGs.tsv")

RIG_with_WUSS_output_dir = 'plots/RIG_WUSS/'

if not os.path.exists(RIG_with_WUSS_output_dir):
    os.makedirs(RIG_with_WUSS_output_dir)

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


def plot_RIG_WUSS(RF, RIG_dict, WUSS_color_dict, encodings, filename='test'):
    """plots RIG - enhanced with WUSS"""

    legend_elements = [Patch(facecolor='r', label='Stem'),
                       Patch(facecolor='g', label='Hairpin Loop'),
                       Patch(facecolor='c', label='Bulge/Interior'),
                       Patch(facecolor='b', label='Others'),
                       Patch(facecolor='m', label='Pseudoknot')]

    # Clean up all pending trials
    plt.clf()

    # Create the figure
    fig, ax = plt.subplots(1, figsize=(20, 10))

    # Plots RIG values
    styles = ['-', '--', '-.', ':']
    for enc, st in zip(encodings, styles[:len(encodings)]):
        ax.plot(RIG_dict[enc][RF], color='k', linewidth=2, ls=st)
        legend_elements.append(
            Line2D([0], [0], ls=st, label=enc, color='k'),
        )

    # Generate the subdivision in colors of the various blocks
    colored_chunks = function_mbr.get_colored_chunks(WUSS_color_dict[RF])

    # For each color block, color under the corresponding area
    for count, colored_chunk in enumerate(colored_chunks):
        plt.fill_between(colored_chunk[0], 1, facecolor=colored_chunk[2], alpha=0.5)
    # plt.fill_between(range(len(RIG_dict[enc][RF])), RIG_dict[enc][RF],
    #               where=list(WUSS_color_dict[RF])=='r',interpolate=True,
    #               facecolor='r', alpha=0.2 )

    ax.set_title(f'-{encodings[0][-2:]}- RIG - {RF}')

    ax.set_xlabel('alignment position (nt)')
    ax.set_ylabel('RIG')
    ax.set_ylim([0, 1])

    # Add the legend
    legend = ax.legend(handles=legend_elements, loc='best', frameon=1)
    frame = legend.get_frame()
    frame.set_facecolor('w')

    plt.tight_layout()
    plt.savefig(filename + ".pdf")
    plt.savefig(filename + ".png", ppi=600)

    plt.close()


for RF_ in RIG_dict[list(RIG_dict.keys())[0]]:
    print(RF_)
    plot_RIG_WUSS(
        RF_,
        RIG_dict=RIG_dict,
        WUSS_color_dict=WUSS_color_dict,
        encodings=['bear90', 'qbear90', 'zbear90'],
        filename=f"{RIG_with_WUSS_output_dir}{RF_}_90"
    )
