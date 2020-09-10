from collections import Counter, defaultdict
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import subprocess

dir_output_RNA_blocks = 'outputs/RNA_Blocks'
dir_output_MBRs = 'outputs/MBRs'
dir_output_blustclust = 'outputs/blustclust'

# funtion to run BlustClust (filter on identity score)
# folder = seq_str_families/sequence/
# identity = threshold expressed by % (ex. 50)
# return a folder with selected sequences
def run_blustClust(basedir_blustclust, folder_seq, identity, name):
    dir_not_similar = os.path.join(dir_output_blustclust, 'not_similar_' + name + '_' + identity)
    if not os.path.exists(dir_not_similar):
        os.makedirs(dir_not_similar)

    blustclust_run = os.path.join(basedir_blustclust, 'blastclust')
    patter_input = os.path.join(folder_seq, '*')
    patter_output = os.path.join(dir_not_similar, '$filename')

    os.system(f'for file in {patter_input}; do filename=$(basename -- $file); echo $filename; {blustclust_run} -i $file -o {patter_output} -p F -S ' + identity + ' ; done')

    print('BlustClust DONE!')
    return dir_not_similar


# function to filter for number of sequences in an RNA family after BlustClust
# folder = folder generated from run_blustClust
# n_seq_family = minimum number of sequences per family
# identity = threshold expressed by % (ex. 50)
def filter_n_seq(folder, n_seq_family, name, identity):
    dir_filter_n_seq = os.path.join(dir_output_blustclust, 'filter_n_seq_' + name + '_' + identity)
    if not os.path.exists(dir_filter_n_seq):
        os.makedirs(dir_filter_n_seq)

    # New folder with only families with more than n_seq_family members
    num_moved = 0
    for fam in sorted(os.listdir(folder)):
        output = subprocess.check_output("wc -l " + os.path.join(folder, fam), shell=True, universal_newlines=True)
        if int(output.strip().split()[0]) >= int(n_seq_family):
            os.system('cp ' + os.path.join(folder, fam) + ' ' + os.path.join(dir_filter_n_seq,  fam))
            num_moved += 1

    # print(f'Moved {num_moved} families with more than {n_seq_family} family members.')

    print('Filter Nseq DONE!')
    return dir_filter_n_seq


# Function to select the bear sequences after BlustClust filtering
def get_bear(folder, folder_bear, name, identity):
    # folder = folder returned by filter_n_seq
    # folder_bear = folder seq_str_families/bear

    dir_bear_filtered = os.path.join(dir_output_blustclust, 'bear_filtered_' + name + '_' + identity)
    if not os.path.exists(dir_bear_filtered):
        os.makedirs(dir_bear_filtered)

    for fam_clean in sorted(os.listdir(folder)):
        seq = []
        with open(os.path.join(folder, fam_clean)) as f:
            for line in f.readlines():
                seq.append(line.split()[0])

        o = open(os.path.join(dir_bear_filtered, fam_clean), "w")
        with open(os.path.join(folder_bear, fam_clean)) as f2:
            line = f2.readline()
            while line:
                if line[0] == ">" and line[1:-1] in seq:
                    o.write(line)
                    line = f2.readline()
                    o.write(line)
                    line = f2.readline()
                    o.write(line)
                    line = f2.readline()
                    o.write(line)
                    line = f2.readline()
                else:
                    line = f2.readline()
        o.close()

    print('Bear sequence research DONE!')
    return dir_bear_filtered


def distributeGaps(gappedReference, ungappedString):
    assert len(gappedReference.replace('-', '')) == len(ungappedString), 'ungapped strings should be equal'
    result = list(ungappedString)
    gaplist = [m.start() for m in re.finditer('-', gappedReference)]

    for gap in gaplist:
        result.insert(gap, '-')
    result = "".join(result)
    return result

# assert(distributeGaps('--abcdef-g-', 'bombasi') == '--bombas-i-')


# Function to Add gap based on RFAM alignment from seed
def add_gap(folder, seed_rfam, name, identity):
    # folder=bear_filtered
    # seed_rfam

    dir_bear_alignment = os.path.join(dir_output_blustclust, 'bear_alignment_' + name + '_' + identity)
    if not os.path.exists(dir_bear_alignment):
        os.makedirs(dir_bear_alignment)

    c = 0
    for fam in sorted(os.listdir(folder)):
        c += 1

        o = open(os.path.join(dir_bear_alignment, fam), "w")
        with open(os.path.join(folder, fam)) as f:
            line = f.readline()
            while line:
                if line[0] == ">":
                    seq_name = line[1:-1]
                    # output = subprocess.check_output("wc -l not_similar/"+fam, shell=True)

                    seq_alignment = subprocess.check_output('zgrep ' + seq_name+' ' + seed_rfam, shell=True, universal_newlines=True)
                    line = f.readline()
                    test = distributeGaps(seq_alignment.split()[1], line[0:-1])
                    line = f.readline()
                    line = f.readline()
                    bear_seq = line
                    gap_pos = []
                    for i, el in enumerate(test):
                        if el == '-':
                            gap_pos.append(i)
                    for i in gap_pos:
                        bear_seq = bear_seq[:i] + "-" + bear_seq[i:]
                    o.write(bear_seq)
                else:
                    line = f.readline()
        o.close()

        # print(fam+"\t"+str(c)+" famiglie su "+str(len(lista_fam_filter2)))

    print('Alignment with gap DONE!')
    print('Bear sequences with gap in {} folder.'.format(dir_bear_alignment))
    return dir_bear_alignment


def decode(bear):
    alph_bear = {'abc': 'a', 'def': 'A', 'ghi=': '=',
          'lmnop': 'l', 'qrstu': 'L', 'vwxyz^': '^',
          '!"#': 'i', '$%&': 'I', '\'()+': '+',
          '234': 'n', '567': 'N', '890>': '>',
          'ABC': 's', 'DEF': 'S', 'GHIJ': '~',
          'KLMN': 'b', 'OPQR': 'B', 'STUVW': '|',
          'YZ~': 'y', '?_|': 'Y', '/\\@': '@',
          '{}[]': '[', ':': ':' , '-': '-'}

    result = ""
    for ch in bear:
        for key in alph_bear:
            if ch in key:
                result += alph_bear[key]
    return result


# Deconding from text (same as decode but from text file)
def decode_from_file(bear, file_name):
    f = open(file_name).readlines()
    alph_bear = [i.split() for i in f]
    alph_bear.append(['-', '-'])
    result = ""
    for ch in bear:
        for group in alph_bear:
            if ch in group[0]:
                result += group[1]
    return result


# Convert from bear to New Alphabet
# folder = bear_alignment(bear_alignment/bear_alignment___62)
def convert_new_bear(folder, name, identity):
    dir_bear_new_alignment = os.path.join(dir_output_blustclust, 'bear_new_alignment_' + name + '_' + identity)
    if not os.path.exists(dir_bear_new_alignment):
        os.makedirs(dir_bear_new_alignment)

    for fam in sorted(os.listdir(folder)):
        o = open(os.path.join(dir_bear_new_alignment, fam), "w")
        with open(os.path.join(folder, fam)) as f:
            line = f.readline()
            while line:
                o.write(decode(line) + "\n")
                line = f.readline()
        o.close()

    print('Decoding BEAR DONE!')


# Convert from bear to New Alphabet (FILE)
# folder = bear_alignment(bear_alignment/bear_alignment___62)
# file_name = alph_mapping.tsv
def convert_new_bear_file(folder, file_name, name, identity):
    dir_bear_new_alignment = os.path.join(dir_output_blustclust, 'bear_new_alignment_' + name + '_' + identity)
    if not os.path.exists(dir_bear_new_alignment):
        os.makedirs(dir_bear_new_alignment)

    for fam in sorted(os.listdir(folder)):
        o = open(os.path.join(dir_bear_new_alignment, fam), "w")
        with open(os.path.join(folder, fam)) as f:
            line = f.readline()
            while line:
                o.write(decode_from_file(line, file_name)+"\n")
                line = f.readline()
        o.close()

    print('Decoding BEAR from File DONE! New alignments are in the {} folder.'.format(dir_bear_new_alignment))


# Create Blocks from bear alignments
def make_blocks(folder, name, identity):
    dir_output = os.path.join(dir_output_RNA_blocks, 'blocks_new_bear_' + name + '_' + identity)
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)

    families = os.listdir(folder)

    for fam in families:
        o = open(dir_output + '/' + fam, "w")
        f = open(folder + fam).readlines()
        zipped = zip(*f)
        v = []
        for col in zipped:
            if '-' not in col:
                c = Counter(col)
                for key in c:
                    if float(c[key])/float(len(col)) > 0.5:
                        v.append(col)
        new = zip(*v)
        for seq in new:
            o.write(''.join(seq))

    print('Created the {} folder.'.format(dir_output))


#Compute expected frequencies
#folder=Blocks/blocks_new_bear_$alph_$id
def expected_frequencies(folder,alph_bear):
    list_=os.listdir(folder)
    for fam in list_:
        f=open(folder+fam)

        line=f.readline()
        while(line):
            c=Counter(line.strip())
            for key in c:
                alph_bear[key].append(c[key])
            line=f.readline()

    for key in alph_bear:
        alph_bear[key]=sum(alph_bear[key])

    #nnew dict with single character frequency
    c=0
    for key in alph_bear:
        c+=alph_bear[key]
    for key in alph_bear:
        alph_bear[key]=float(alph_bear[key])/float(c)
    v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
    fr_expected=pd.DataFrame(columns=v_bear, index=v_bear)

    for index, row in fr_expected.iterrows():
        for col in v_bear:
            if index==colonna:
                fr_expected.ix[index, col] = alph_bear[index]*2
            else:
                fr_expected.ix[index, col] = 2*alph_bear[index]*alph_bear[col]

    fr_expected.to_csv('expected_frequencies.tsv', sep="\t")
    # print('Expected_frequencies DONE!')
    return fr_expected


# observed_substitution(f_ij)
# folder=Blocks/blocks_new_bear_$alph_$id
# v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
def observed_substitution(dir_output_MBRs_name_id, folder, v_bear, name, identity):
    substitution = pd.DataFrame(1.0, columns=v_bear, index=v_bear)

    list_ = sorted(os.listdir(folder))
    for fam in list_:
        print(fam)
        with open(os.path.join(folder, fam)) as f:
            f = [x.strip() for x in f.readlines()]
            v1, v2 = np.triu_indices(len(f))

            for k in range(0, len(v1)):
                if v1[k] != v2[k]:
                    for i, el in enumerate(f[v1[k]]):
                        tmp = f[v2[k]][i]
                        substitution.at[el, tmp] += 1.0

    substitution = substitution + substitution.T - 1.0
    substitution.to_csv(os.path.join(dir_output_MBRs_name_id, 'substitution_' + name + '_' + identity + '.tsv'), sep="\t")

    print('Substitution matrix DONE!')
    return substitution


# substitution=dataframe returned by observed_substitution
# v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
def make_q(dir_output_MBRs_name_id, substitution, v_bear, name, identity):
    number_couple = 0
    for j in range(0, len(v_bear)):
        number_couple += substitution.iloc[j, j:].sum()
    # print(number_couple)

    q_ij = substitution.divide(number_couple)
    q_ij.to_csv(os.path.join(dir_output_MBRs_name_id, 'q_ij_' + name + '_' + identity + '.tsv'), sep="\t")

    print('q_ij DONE!')
    return q_ij


# q_ij returned by make_q
# v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
def make_p(q_ij, v_bear):
    p_i = dict.fromkeys(v_bear, 0)
    for char in v_bear:
        same = q_ij[char][char]
        others = (sum(q_ij[char]) - same) / 2
        p_i[char] += same + others

    print('p_i dict DONE!')
    return p_i


# p_i returned by make_p
# v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
def make_e(dir_output_MBRs_name_id, p_i, v_bear, name, identity):
    e_ij = pd.DataFrame(1.0, columns=v_bear, index=v_bear)

    for char in v_bear:
        for char2 in v_bear:
            if char == char2:
                e_ij.loc[char, char2] = p_i[char] * p_i[char2]
            else:
                e_ij.loc[char, char2] = 2 * p_i[char] * p_i[char2]

    e_ij.to_csv(os.path.join(dir_output_MBRs_name_id, 'E_ij_' + name + '_' + identity + '.tsv'), sep='\t')

    print('e_ij DONE!')
    return e_ij


# freq_observed=dataframe of observed frequencies (q_ij)
# fr_expected=dataframe of expected frequencies (e_ij)
def make_matrix(dir_output_MBRs_name_id, freq_observed, fr_expect, name, identity):
    # Odds ratio matrix
    odds_ratio_matrix = freq_observed.divide(fr_expect)

    odds_ratio_matrix.to_csv(os.path.join(dir_output_MBRs_name_id, 'odds_ratio_matrix_' + name + '_' + identity + '.tsv'), sep="\t")
    print('Odds ratio matrix DONE!')

    # Score Matrix
    mbr_new = odds_ratio_matrix.applymap(np.log2)
    mbr_new.to_csv(os.path.join(dir_output_MBRs_name_id, 'MBR_' + name + '_' + identity + '.tsv'), sep="\t")
    print('MBRs DONE!')
    return mbr_new


def Expected_score(s_ij, p_i):
    E = 0
    for i, char in enumerate(s_ij.columns.values):
        for i2, char2 in enumerate(s_ij.columns.values):
            if i2 >= i:
                # print(char, char2)
                E += s_ij.loc[char, char2] * p_i[char] * p_i[char2]
    # print(E)
    return E


def entropy(q_ij, s_ij):
    H = 0
    length = len(q_ij.columns.value_counts())
    for j in range(0, length):
        H += (q_ij.iloc[j, j:] * s_ij.iloc[j, j:]).sum()

    # print(H)
    return H


def make_heatmap(dir_output_MBRs_name_id, S_ij, name, identity, encoding_size):
    sns.set(font_scale=2.0)
    plt.figure(figsize=(encoding_size, encoding_size))
    sns.heatmap(S_ij, xticklabels=1, yticklabels=1, cmap="YlGnBu")
    plt.savefig(os.path.join(dir_output_MBRs_name_id, 'matrix_' + name + '_' + identity + '.pdf'))
    plt.close()


# Function that start from sequence/structure of all RFAM families and create a new MBRs from alphabet file
# folder = folder with RFAM RNA sequences
# folder_bear = folder with RFAM RNA structure in bear
# RFAM_seed_file = Rfam seed downloaded from RFAM database
# id_blustClust = Sequence identity for filter
# filter_nSeq = threshold on number of sequences in a family after filtering
# file_alph = alphabet file with bear mapping
# file_info = output file with information about the MBRs

def BlustClust_filter_alignment(basedir_blustclust, folder_seq, folder_bear, RFAM_seed_file_gz, id_blustClust, filter_nSeq, file_alph):
    name = os.path.basename(file_alph).split('.')[0]
    dir_not_similar = run_blustClust(basedir_blustclust, folder_seq, id_blustClust, name)
    dir_filter_n_seq = filter_n_seq(dir_not_similar, filter_nSeq, name, id_blustClust)
    dir_bear_filtered = get_bear(dir_filter_n_seq, folder_bear, name, id_blustClust)
    dir_bear_alignment = add_gap(dir_bear_filtered, RFAM_seed_file_gz, name, id_blustClust)

    # convert_new_bear_file(os.path.join('bear_alignment___'+id_blustClust), file_alph, name, id_blustClust)
    convert_new_bear_file(dir_bear_alignment, file_alph, name, id_blustClust)


def Make_MBR_from_blocks(blocks_folder, id_blustClust, file_alph, file_info):
    name = os.path.basename(file_alph).split('.')[0]

    dir_output_MBRs_name_id = os.path.join(dir_output_MBRs, name + '_' + id_blustClust)
    if not os.path.exists(dir_output_MBRs_name_id):
        os.makedirs(dir_output_MBRs_name_id)

    with open(file_alph) as f:
        v_char = [x.split()[1] for x in f.readlines()]

    sost = observed_substitution(dir_output_MBRs_name_id, blocks_folder, v_char, name, id_blustClust)
    q_ij = make_q(dir_output_MBRs_name_id, sost, v_char, name, id_blustClust)  # freq_observed

    p_i = make_p(q_ij, v_char)
    e_ij = make_e(dir_output_MBRs_name_id, p_i, v_char, name, id_blustClust)   # fr_expect

    S_ij = make_matrix(dir_output_MBRs_name_id, q_ij, e_ij, name, id_blustClust)

    E = Expected_score(S_ij, p_i)
    H = entropy(q_ij, S_ij)

    with open(os.path.join(dir_output_MBRs_name_id, file_info), "w") as fw:
        fw.write('Expected_score:\t' + str(E) + '\nEntropy:\t' + str(H))

    make_heatmap(dir_output_MBRs_name_id, S_ij, name, id_blustClust, len(v_char))


def get_colored_chunks(color_string, base1=False):
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
            rangeY = [1, 1]

        if base1:
            rangeX = rangeX + 1

        colored_chunks.append([rangeX, rangeY, my_color])

    return colored_chunks


def load_and_color_WUSS_dictionary(WUSS_path):
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
        color_str = ""
        for c in WUSS:
            if c in stems:
                color_str += "r"
            elif c in hairpinloops:
                color_str += "g"
            elif c in bulgeInterior:
                color_str += "c"
            elif c in gaps:
                color_str += "w"
            elif c.isalpha():
                color_str += "m"
            else:
                color_str += "b"

        WUSS_color_dict[RF] = color_str

    return WUSS_color_dict


def load_rig_values(path_encoding_list):
    RIG_dict = defaultdict(dict)

    for path, encoding in path_encoding_list:
        with open(path) as f:
            for line in f:
                data = line.split()
                RIG_dict[encoding][data[0]] = [float(x) for x in data[1:]]

    return RIG_dict
