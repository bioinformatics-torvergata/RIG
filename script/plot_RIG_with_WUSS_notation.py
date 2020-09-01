# RIG plots with WUSS notation from secondary structure consensus.

# Background colors are mapped to secondary structure elements of the consensus structure provided by the covariance
# model.

from collections import defaultdict

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import seaborn as sns
sns.set(context='paper', style='white', palette='deep', font='serif', font_scale=2, color_codes=True, rc=None)


# Paths and directories
WUSS_path = "SS_cons/SS_cons_WUSS.tsv"
RIG_dir = "RIG_nogaps/"  # gaps removed from the consensus

bear90_path = RIG_dir + "bear_90_RIGs.tsv"
bear50_path = RIG_dir + "bear_50_RIGs.tsv"
qbear90_path = RIG_dir + "qbear_90_RIGs.tsv"
qbear50_path = RIG_dir + "qbear_50_RIGs.tsv"
zbear90_path = RIG_dir + "zbear_90_RIGs.tsv"
zbear50_path = RIG_dir + "zbear_90_RIGs.tsv"

RIG_with_WUSS_output_dir = 'RIG_WUSS/'

# Symbols
stems = "><}{][)("
hairpinloops = "_"
bulgeInterior = "-"
gaps = "."

# Read WUSS dictionary
WUSS_dict = dict()
with open(WUSS_path) as f:
    for line in f:
        # RF\tWUSS
        RF, WUSS = line.split('\t')

        # Remove gaps from the consensus
        WUSS_dict[RF] = WUSS.strip().replace(gaps, "")


# Colors WUSS dictionary
WUSS_color_dict = dict()

for RF, WUSS in WUSS_dict.items():
    colorStr = ""
    for c in WUSS:
        if c in stems:
            colorStr += "r"
        elif c in hairpinloops:
            colorStr += "g"
        elif c in bulgeInterior:
            colorStr += "c"
        elif c in gaps:
            colorStr += "w"
        elif c.isalpha():
            colorStr += "m"
        else:
            colorStr += "b"

    WUSS_color_dict[RF] = colorStr


RIG_dict = defaultdict(dict)


# Load RIG values
def load_rig_values(path, encoding):
    with open(path) as f:
        for line in f:
            data = line.split()
            RF = data[0]
            rig = [float(x) for x in data[1:]]
            RIG_dict[encoding][RF] = rig


load_rig_values(bear90_path, 'bear90')
load_rig_values(bear50_path, 'bear50')
load_rig_values(qbear90_path, 'qbear90')
load_rig_values(qbear50_path, 'qbear50')
load_rig_values(zbear90_path, 'zbear90')
load_rig_values(zbear50_path, 'zbear50')


def get_colored_chunks(RIG_values, color_string):
    """assigns color values to the SSE"""

    i = 0
    colored_chunks = []
    while i < len(color_string):
        # Starting color
        my_color = color_string[i]
        start_idx = i
        while (i < len(color_string)) and (color_string[i] == my_color):
            i += 1

        end_idx = i

        # rangeX = np.arange(start_idx, end_idx, 1)
        rangeX = np.arange(start_idx - 0.1, end_idx + 0.1)

        # rangeY = RIG_values[start_idx:end_idx]
        rangeY = 1 * np.arange(start_idx, end_idx, 1)

        # If it is a single point then draw a narrow band at point height
        if end_idx - start_idx == 1:
            # rangeX = [start_idx-0.2, start_idx+0.2]
            # rangeY = 2*[RIG_values[start_idx]]
            rangeY = [1, 1]

        colored_chunks.append([rangeX, rangeY, my_color])

    return colored_chunks


def plot_RIG_WUSS(RF, RIG_dict, WUSS_color_dict, filename='test', encodings=None):
    """plots RIG - enhanced with WUSS"""

    if encodings is None:
        encodings = ['bear90']

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
    colored_chunks = get_colored_chunks(RIG_dict[enc][RF], WUSS_color_dict[RF])

    # For each color block, color under the corresponding area
    for count, colored_chunk in enumerate(colored_chunks):
        plt.fill_between(colored_chunk[0], 1, facecolor=colored_chunk[2], alpha=0.5)
    # plt.fill_between(range(len(RIG_dict[enc][RF])), RIG_dict[enc][RF],
    #               where=list(WUSS_color_dict[RF])=='r',interpolate=True,
    #               facecolor='r', alpha=0.2 )

    ax.set_title(f'-{enc[-2:]}- RIG - {RF}')

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


for RF_ in RIG_dict['bear90']:
    print(RF_)
    plot_RIG_WUSS(RF_,
                  RIG_dict=RIG_dict,
                  WUSS_color_dict=WUSS_color_dict,
                  encodings=['bear90', 'qbear90', 'zbear90'],
                  filename=f"{RIG_with_WUSS_output_dir}{RF_}_90"
                  )
