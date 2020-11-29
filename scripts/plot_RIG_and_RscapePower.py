import pandas as pd
import numpy as np
import os
import re
import sys

import function_mbr

import matplotlib.pylab as plt
import seaborn as sns

sns.set()
sns.set_context(rc={"font.size": 0, "axes.titlesize": 34, "axes.labelsize": 25, "legend.fontsize": 20})


WUSS_path = sys.argv[1]
dir_input_rscape = sys.argv[2]

dir_output_RIG_RscapePower_plots = 'outputs/plots/RIG_RscapePower/'
RIG_and_RscapePower_tsv_path = os.path.join(dir_output_RIG_RscapePower_plots, 'RIG_and_Rscape.tsv')

if not os.path.exists(dir_output_RIG_RscapePower_plots):
    os.makedirs(dir_output_RIG_RscapePower_plots)


RIG_dir = "outputs/RIGs/"

bear90_path = os.path.join(RIG_dir, "bear_90_RIGs.tsv")
bear50_path = os.path.join(RIG_dir, "bear_50_RIGs.tsv")
qbear90_path = os.path.join(RIG_dir, "qbear_90_RIGs.tsv")
qbear50_path = os.path.join(RIG_dir, "qbear_50_RIGs.tsv")
zbear90_path = os.path.join(RIG_dir, "zbear_90_RIGs.tsv")
zbear50_path = os.path.join(RIG_dir, "zbear_50_RIGs.tsv")

# Load RIG values
RIG_dict = function_mbr.load_rig_values([
    [bear90_path, 'bear90'],
    # [bear50_path, 'bear50'],
    [qbear90_path, 'qbear90'],
    # [qbear50_path, 'qbear50'],
    [zbear90_path, 'zbear90'],
    # [zbear50_path, 'zbear50'],
])

one_read_encoding = list(RIG_dict.keys())[0]

# Parsing to extract the R-scape power from R-scape output.
RFXXX_noInfo_list = []
RFXXX_missingRIGs_list = []
for xxx_txt in sorted([x for x in os.listdir(dir_input_rscape) if x.endswith('.txt')]):
    path_xxx_txt = os.path.join(dir_input_rscape, xxx_txt)
    RFXXX = path_xxx_txt.split('/')[-1].split('.')[0]

    if RFXXX not in RIG_dict[one_read_encoding]:
        print('Missing RIG scores for', RFXXX, 'family')
        RFXXX_missingRIGs_list.append(RFXXX)
        continue

    with open(path_xxx_txt) as f:
        file_splitted_list = f.read().split('#----------------------------------------------------------------')
        if len(file_splitted_list) != 4:
            # print(file_splitted_list)
            RFXXX_noInfo_list.append(RFXXX)
            continue

        RIG_dict['rscape'][RFXXX] = [-1] * len(RIG_dict[one_read_encoding][RFXXX])
        for line in file_splitted_list[3].strip('\n').split('\n'):
            if not line.startswith('#'):
                # print(line)
                left_pos, right_pos, substitutions, power = re.findall(
                    "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line
                )
                # print(line, '-->', left_pos, right_pos, substitutions, power)

                for xxx_pos in [left_pos, right_pos]:
                    RIG_dict['rscape'][RFXXX][int(xxx_pos) - 1] = float(power)

print('Num. Rfam families with R-scape power analysis information:', len(RIG_dict['rscape']))
# print(len(RFXXX_noInfo_list) + len(RFXXX_missingRIGs_list))

# Remove gaps
if True:
    with open(WUSS_path) as f:
        for line in f:
            RF, WUSS = line.strip('\n').split('\t')

            if RF in RIG_dict[one_read_encoding]:
                for encoding, rf_to_rigs_dict in RIG_dict.items():
                    if encoding != 'rscape' or RF in RIG_dict['rscape']:
                        RIG_dict[encoding][RF] = function_mbr.removeGapsFromRIG(WUSS, rf_to_rigs_dict[RF])
                    # else:
                    #    print('Missing Rscape values for', RF, 'family')
            # else:
            #    print('Missing RIG scores for', RF, 'family')

# For each encoding there is a table of RIG scores: each row is a Rfam family, each column is a position in the
# alignment of that family.


max_len = 0
for rf_to_rigs_dict in RIG_dict.values():
    for rig_list in rf_to_rigs_dict.values():
        if len(rig_list) > max_len:
            max_len = len(rig_list)
# print('max_len', max_len)

RIGS_dfs = {}
for encoding, rf_to_rigs_dict in RIG_dict.items():
    for rf, rig_list in rf_to_rigs_dict.items():
        if len(rig_list) < max_len:
            RIG_dict[encoding][rf] = rig_list + [np.nan] * (max_len - len(rig_list))

    RIGS_dfs[encoding] = pd.DataFrame.from_dict(RIG_dict[encoding], orient='index')

# Parsing to have all the information in one file.
if not os.path.exists(RIG_and_RscapePower_tsv_path):
    with open(RIG_and_RscapePower_tsv_path, 'w') as fw:
        fw.write('\t'.join(['encoding', 'family_code', 'position', 'score']) + '\n')

        to_append = ''
        for xxx_xx, dataframe in RIGS_dfs.items():
            if xxx_xx != 'rscape':
                c = 0
                for RFXXX in dataframe.index:
                    if RFXXX in RIG_dict['rscape']:
                        for i, x in enumerate(dataframe.loc[RFXXX].dropna()):
                            pos = i
                            power = RIG_dict['rscape'][RFXXX][pos]
                            fw.write('\t'.join([str(x) for x in [xxx_xx, RFXXX, pos, x]]) + '\n')

                            # Only for one encoding (it is the same for all the others)
                            if xxx_xx == one_read_encoding:
                                to_append += '\t'.join([str(x) for x in ['rscape', RFXXX, pos, power]]) + '\n'

                # print(xxx_xx, dataframe.min().min(), dataframe.max().max(), c)

        fw.write(to_append)

# Read the created file with RIG score and R-scape power values.
rig_and_rscape_df = pd.read_csv(RIG_and_RscapePower_tsv_path, sep='\t')

# NOTE: to execute if you don't want to visualize the RIG scores calculated using RNA Blocks obtained by removing
# redundant primary sequences up to 50% of similarity for each Rfam seed alignment.
rig_and_rscape_df = rig_and_rscape_df[rig_and_rscape_df['encoding'].str.endswith('_50') == False]

# Generate plots
y_min, y_max = -0.01, 1.01

for RFXXX in rig_and_rscape_df['family_code'].unique():
    path_image = os.path.join(dir_output_RIG_RscapePower_plots, RFXXX + '.rscape.png')
    print(path_image)
    if not os.path.exists(path_image):
        # print(RFXXX, rig_and_rscape_df[(rig_and_rscape_df['family_code'] == RFXXX)]['position'].max())

        try:
            ax = sns.pointplot(
                data=rig_and_rscape_df[(rig_and_rscape_df['family_code'] == RFXXX)], x='position', y='score',
                hue='encoding', fit_reg=False,
                markers=["o", "o", "o", "x"],
                linestyles=["-", "-", "-", "--"]
            )

            ax.set_title(RFXXX)

            ax.set_xlim(-0.5, rig_and_rscape_df[(rig_and_rscape_df['family_code'] == RFXXX)]['position'].max() + 0.5)
            ax.set_ylim(y_min, y_max)

            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, size=20)
            # ERROR:ax.set_yticklabels(['{:.2f}'.format(float(t.get_text().replace('−', '-'))) for t in ax.get_yticklabels()], size = 25)
            ax.set_yticklabels(ax.get_yticklabels(), size=25)

            ax.set_yticks([y_min, 0, 0.25, 0.5, 0.75, 1, y_max])
            ax.set_yticklabels(['', '0', '0.25', '0.5', '0.75', '1', ''], size=25)  # set the labels

            fig = ax.get_figure()
            fig.set_size_inches(
                rig_and_rscape_df[(rig_and_rscape_df['family_code'] == RFXXX)]['position'].max() / 3, 10
            )

            # plt.setp(ax.get_legend().get_texts(), fontsize='22') # for legend text
            plt.setp(ax.get_legend().get_title(), fontsize='24')  # for legend title

            fig.tight_layout()

            for x in rig_and_rscape_df[
                (rig_and_rscape_df['family_code'] == RFXXX) & (rig_and_rscape_df['encoding'] == 'rscape')][['position', 'score']
            ].iterrows():
                position, rscape_power = x[1]
                position += 1
                if rscape_power < 0:
                    # ax.axvspan(xmin=position-2, xmax=position-1, ymin=-1, ymax=1, facecolor='b', alpha=0.5)
                    ax.axvspan(xmin=position - 1.5, xmax=position - 0.5, ymin=-1, ymax=1, facecolor='g', alpha=0.5)
            fig.savefig(path_image)
            fig.clf()
        except:
            print('Problem with ' + RFXXX)
