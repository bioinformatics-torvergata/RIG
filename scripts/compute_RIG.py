import os
import sys

import numpy as np
from scipy.stats import entropy

import pickle

dir_output_RIGs = 'outputs/RIGs'


def rescaleColumn(col, **kwargs):
    """col is an array of floats
    it outputs a rescaled column
    Needed as a reminder that we can try different normalization (PSSM to PPM transition for IC)
    """
    if kwargs.get('rescale', 0):
        if kwargs.get('rescale') == "delogit":
            """exp(logodds) -> o/(1+o) -> /sum()"""
            odds = np.exp(col)
            delogp = odds / (1 + odds)
            # somma fissa diversa da 1, credo dipenda dalla matrice di sostituzione
            return delogp / delogp.sum()
        if kwargs.get('rescale') == "delogit2":
            """exp(logodds) -> /sum() qualitative arguments"""
            odds = np.exp(col)

            # somma fissa diversa da 1, credo dipenda dalla matrice di sostituzione
            return odds / odds.sum()
        print(
            f"[WARNING] This rescale method ('{kwargs.get('rescale')}') is not yet implemented! Using standard rescaling")
        tmp = np.array(col)
        M = tmp.max()
        m = tmp.min()
        return (tmp - m) / (M - m)
    else:  ##base case, rescale to 1, no transform
        tmp = np.array(col)
        M = tmp.max()
        m = tmp.min()
        return (tmp - m) / (M - m)


def IC(col):
    """the entropy method of scipy stats automatically normalize the distribution in case it does not sum to 1
    returns (maxE-E)/maxE : the fraction of maximum entropy in that PSSM aka the information gain normalized from
    0 to 1 (i.e. it is 0 if there is max entropy, all characters, and 1 if only one character is present, in
    this way it is comparable between alphabets of different sizes)"""
    maxE = np.log2(col.shape)
    E = entropy(col, qk=None, base=2)
    return (maxE - E) / maxE


def computePSSMIC(PSSM, **kwargs):
    """applies rescaling (if different from standard) and scipy.stats.entropy to every column of PSSM"""
    ic_arr = []
    for col in PSSM:
        col = rescaleColumn(col, rescale=kwargs.get('rescale', None))
        ic_arr.append(IC(col))

    return np.array(ic_arr)


def computeAllIC(PSSMs, **kwargs):
    ICarrays = []
    tot = len(PSSMs)
    for i, PSSM in enumerate(PSSMs):
        print(f"doing PSSM {i + 1} of {tot}")
        ICarrays.append(computePSSMIC(PSSM, **kwargs))
    return np.array(ICarrays)


xxx_pssm_path = sys.argv[1]

alphaexts=[]
path_pssm_list = [xxx_pssm_path]

for ps in path_pssm_list:
    filename = ps.split("/")[-1].split(".")[0]
    alphaexts.append("_".join(filename.split("_")[-2:]))

PSSMs = []
rflists = []
for path_pssm in path_pssm_list:
    with open(path_pssm, 'rb') as f1:
        tmp = pickle.load(f1)

    if isinstance(tmp, dict):
        klist = list(tmp.keys())
        PSSMs.append([tmp[rf] for rf in klist])
        rflists.append(klist)
    else:
        print(f"I should not be here, are you sure the PSSMs is a huge dict and not just a {type(tmp)}?")
        PSSMs.append(tmp)
        # load associated rfam list?? better build everything as dictionary

ICarrays = []
rescale = 'delogit2'
for PSSMset in PSSMs: #diversi file da testare
    #must transpose to cycle over columns
    PSSMsetT = [profile.T for profile in PSSMset]
    ICarrays.append(computeAllIC(PSSMsetT, rescale=rescale)) #this will have shape (#RF, )


if not os.path.exists(dir_output_RIGs):
    os.makedirs(dir_output_RIGs)

###write file with RIG,
## make a DF for each file
# the table shall be written

for fileno in range(len(path_pssm_list)):
    myenc = alphaexts[fileno]
    filename = f"{myenc}_RIGs.tsv"
    with open(os.path.join(dir_output_RIGs, filename), 'w') as fw:
        for rf, rigs in zip(rflists[fileno], ICarrays[fileno]):
            rigrow = "\t".join([str(el[0]) for el in rigs])
            row = f"{rf}\t{rigrow}\n"
            fw.write(row)
