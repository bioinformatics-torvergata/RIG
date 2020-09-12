import pandas as pd
import os
import re

import matplotlib.pylab as plt
import seaborn as sns

sns.set()
sns.set_context(rc={"font.size": 0, "axes.titlesize": 34, "axes.labelsize": 25, "legend.fontsize": 20})

dir_output_rscape = 'data/Rfam_stockholm_rscapes'

RIG_sheet_path = 'data/RIG/withgaps/All_RIGs.filled_columns.xlsx'

dir_RIG_and_RscapePower = 'scripts/RIG_RscapePower/'
RIG_and_RscapePower_tsv_path = os.path.join(dir_RIG_and_RscapePower, 'RIG_and_Rscape.tsv')

dir_output_RIG_RscapePower_plots = 'plots/RIG_RscapePower/'

if not os.path.exists(dir_output_RIG_RscapePower_plots):
    os.makedirs(dir_output_RIG_RscapePower_plots)


# Read RIG scores
RIGS_dfs = pd.read_excel(RIG_sheet_path, sheet_name=None, header=None)


# For each encoding there is a table of RIG scores: each row is a Rfam family, each column is a position in the
# alignment of that family.
for xxx_xx, dataframe in RIGS_dfs.items():
    RIGS_dfs[xxx_xx] = dataframe.set_index(0)


# Parsing to extract the R-scape power from R-scape output.
RFXXX_to_pos_to_power_dict = {}

c_list = []
for xxx_txt in sorted(os.listdir(dir_output_rscape)):
    path_xxx_txt = os.path.join(dir_output_rscape, xxx_txt)
    RFXXX = path_xxx_txt.split('/')[-1].split('.')[0]
    # print(path_xxx_txt, RFXXX)

    with open(path_xxx_txt) as f:
        file_splitted_list = f.read().split('#----------------------------------------------------------------')
        if len(file_splitted_list) != 4:
            # print(file_splitted_list)
            c_list.append(RFXXX)
            continue

        RFXXX_to_pos_to_power_dict[RFXXX] = {}

        for line in file_splitted_list[3].strip('\n').split('\n'):
            if not line.startswith('#'):
                # print(line)
                left_pos, right_pos, substitutions, power = re.findall(
                    "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line
                )
                # print(line, '-->', left_pos, right_pos, substitutions, power)

                for xxx_pos in [left_pos, right_pos]:
                    RFXXX_to_pos_to_power_dict[RFXXX][int(xxx_pos)] = power

print('Num. Rfam families with R-scape power analysis information:', len(RFXXX_to_pos_to_power_dict))


# Parsing to have all the information in one file.
if not os.path.exists(RIG_and_RscapePower_tsv_path):
    if not os.path.exists(dir_RIG_and_RscapePower):
        os.makedirs(dir_RIG_and_RscapePower)

    with open(RIG_and_RscapePower_tsv_path, 'w') as fw:
        fw.write('\t'.join(['encoding', 'family_code', 'position', 'score']) + '\n')

        to_append = ''
        for xxx_xx, dataframe in RIGS_dfs.items():
            c = 0
            for RFXXX in dataframe.index:
                if RFXXX in RFXXX_to_pos_to_power_dict:
                    for i, x in enumerate(dataframe.at[RFXXX].dropna()):
                        pos = i + 1
                        power = -1
                        if pos in RFXXX_to_pos_to_power_dict[RFXXX]:
                            power = RFXXX_to_pos_to_power_dict[RFXXX][pos]
                        fw.write('\t'.join([str(x) for x in [xxx_xx, RFXXX, pos, x]]) + '\n')

                        # Only for one encoding (it is the same for all the others)
                        if xxx_xx == 'bear_50':
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

            ax.set_xlim(-0.5, rig_and_rscape_df[(rig_and_rscape_df['family_code'] == RFXXX)]['position'].max() - 0.5)
            ax.set_ylim(y_min, y_max)

            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, size=20)
            # ERROR:ax.set_yticklabels(['{:.2f}'.format(float(t.get_text().replace('−', '-'))) for t in ax.get_yticklabels()], size = 25)
            ax.set_yticklabels(ax.get_yticklabels(), size=25)

            ax.set_yticks([y_min, 0, 0.25, 0.5, 0.75, 1, y_max])
            ax.set_yticklabels(['', '0', '0.25', '0.5', '0.75', '1', ''], size=25)  # set the labels

            fig = ax.get_figure()
            fig.set_size_inches(rig_and_rscape_df[(rig_and_rscape_df['family_code'] == RFXXX)]['position'].max() / 3,
                                10)

            # plt.setp(ax.get_legend().get_texts(), fontsize='22') # for legend text
            plt.setp(ax.get_legend().get_title(), fontsize='24')  # for legend title

            fig.tight_layout()

            for x in rig_and_rscape_df[
                (rig_and_rscape_df['family_code'] == RFXXX) & (rig_and_rscape_df['encoding'] == 'rscape')][
                ['position', 'score']].iterrows():
                position, rscape_power = x[1]
                # print(position, rscape_power)
                if rscape_power < 0:
                    # ax.axvspan(xmin=position-2, xmax=position-1, ymin=-1, ymax=1, facecolor='b', alpha=0.5)
                    ax.axvspan(xmin=position - 1.5, xmax=position - 0.5, ymin=-1, ymax=1, facecolor='g', alpha=0.5)
            fig.savefig(path_image)
            fig.clf();
        except:
            print('Problem with ' + RFXXX)
