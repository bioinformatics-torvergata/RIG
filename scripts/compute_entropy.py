import os
from collections import defaultdict
from collections import Counter
import gzip
import numpy as np
import sys
import pickle

# Paths and directories
dir_output_entropy = 'outputs/entropy/'

if not os.path.exists(dir_output_entropy):
    os.makedirs(dir_output_entropy)

WUSS_path = sys.argv[1]

path_gapfam_pickle_gz = sys.argv[2]
with gzip.open(path_gapfam_pickle_gz, 'rb') as afile:
    gapped_fam_dict = pickle.load(afile)

# RF from consensus list and request primary sequence alignments from RFAM
def fromTextToAlign(seqId_to_info_to_str_dict):
    al = []
    for seq_id, info_to_str_dict in seqId_to_info_to_str_dict.items():
        al.append(info_to_str_dict['sequence'])
    return al


# Read WUSS dictionary
families = defaultdict(dict)

with open(WUSS_path) as f:
    for line in f:
        # RF\tWUSS
        rf_, cons = line.strip('\n').split('\t')
        families[rf_]['SS_cons'] = cons.strip()

        families[rf_]['align'] = fromTextToAlign(gapped_fam_dict[rf_])


# 2 remove gaps by looking at the consensus
def removeGapsFromAlign(sscons, align, gapCharacter="."):
    """
    removes positions from the alignment that match gaps in the sscons
    input:
        sscons in WUSS format
        alignment as string array
    output:
        modified alignment
    """

    assert len(sscons) == len(align[0]), "different lengths!\n sscons w gaps:{}, align:{}".format(
        (len(sscons)), len(align[0])
    )

    cleanAlign = ["".join([a for ss, a in zip(sscons, seq) if ss != gapCharacter]) for seq in align]
    if (len(sscons) - sscons.count(gapCharacter)) != len(cleanAlign[0]):
        raise ValueError('gaps removal error!')

    return cleanAlign


# print('test clean align')
# print(removeGapsFromAlign("aaaa.....aaaa", ["AAAABBBBBAAAA","AAAABBBBBAAAA","BBBBAAAAABBBB"]))

for rf_ in families:
    # print(rf_)
    errors = []
    try:
        families.get(rf_)['gapless'] = removeGapsFromAlign(
            families.get(rf_).get('SS_cons'), families.get(rf_).get('align')
        )
    except AssertionError as ae:
        print(rf_, ae)
        errors.append(rf_)


def alignmentToShannon(align):
    """
    given an array of strings, compute the per position shannon entropy
    input:
        - array of strings (same length)
    output:
        - array of shannon entropies
    """

    if align is None:
        return []

    shannon = []
    for pos in range(len(align[0])):
        column_string = "".join([seq[pos] for seq in align])
        H = seq_RIG_H(column_string)
        shannon.append(H)

    return shannon


def seq_RIG_H(s, char_amount=4):
    """log2 entropy of a string (excluding gaps)
    char_amount is used to get the normalized RIG (4 if standard nt)
    """

    gap_proportion = s.count("-") / len(s)
    s = s.replace("-", "")
    probabilities = [n_x / len(s) for x, n_x in Counter(s).items()]
    e_x = [-p_x * np.log2(p_x) for p_x in probabilities]

    # Gaps proportionally reduce the final score
    # (as per https://onlinelibrary.wiley.com/doi/full/10.1002/prot.10146)
    return (np.log2(char_amount) - sum(e_x)) / np.log2(char_amount) * (1 - gap_proportion)


# print('test seq_RIG_H')
# print(seq_RIG_H('aaaaaaaab'))

# print('test alignmentToShannon')
# print(alignmentToShannon([
#    'GAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUC---AACUCGAA-UAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCC---',
#    'GAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUC---AACUCGAA-UAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCC---',
#    'GAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUC---AACUCGAA-UAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCCGCU',
#    'AGAACGGAGCGGUUUCUCGUUUAACCCUUGAAGACACCGCCCGUUCAGAGGGUAUCUCUCGAACCCGAAAUAACUAAAGCCAACGUGAACUUUUGCGGACCUC--UGGUCCGCU']))


num_families_div_10 = len(families) // 10

# Compute Entropy
for num_fam, rf_ in enumerate(families):
    if num_fam % num_families_div_10 == 0:
        print('{:.2f}% ({}/{})'.format((num_fam / len(families)) * 100.0, num_fam, len(families)))

    errorsEntropy = []
    try:
        families.get(rf_)['entropy'] = alignmentToShannon(families.get(rf_).get('gapless', None))
    except AssertionError:
        errorsEntropy.append(rf_)


def printEntropy(family_dict, fname="test"):
    """given a dict d[rf]['entropy'] print on file"""

    with open(fname, 'w') as fw:
        for rf_ in family_dict:
            line = rf_ + "\t" + "\t".join([str(x) for x in family_dict[rf_]['entropy']]) + "\n"
            fw.write(line)

    print("File '{}' created.".format(fname))


printEntropy(families, os.path.join(dir_output_entropy, "entropy.tsv"))
