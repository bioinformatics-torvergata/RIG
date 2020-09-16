import os
import gzip
import subprocess
import re
import pickle
import sys

import function_mbr

path_rfam_seed_gz = sys.argv[1]
dir_base = sys.argv[2]
path_RNAfold = sys.argv[3]
path_BEAR = sys.argv[4]

dir_SS_cons = os.path.join(dir_base, 'SS_cons')
dir_gapped_fam = os.path.join(dir_base, 'gapped_fam')
dir_rfam_stockholm = os.path.join(dir_base, 'stockholm')
dir_rfam_constraint = os.path.join(dir_base, 'constraint')
dir_rfam_folded = os.path.join(dir_base, 'folded')
dir_rfam_bear = os.path.join(dir_base, 'bear')
dir_rfam_sequences = os.path.join(dir_base, 'sequences')

for dir_xxx in [dir_SS_cons, dir_gapped_fam, dir_rfam_stockholm, dir_rfam_constraint, dir_rfam_folded, dir_rfam_folded,
                dir_rfam_bear, dir_rfam_sequences]:
    if not os.path.exists(dir_xxx):
        os.makedirs(dir_xxx)

path_SS_cons_WUSS = os.path.join(dir_SS_cons, 'SS_cons_WUSS.tsv')
path_GC_RF = os.path.join(dir_SS_cons, 'GC_RF.tsv')
path_gapped_fam_dict_pickle = os.path.join(dir_gapped_fam, 'gapped_fam_dict.pickle.gz')

# Detect duplicated IDs
dupId_to_freq_dict = {}
dupId_to_suffix_dict = {}
ids_set = set()
with gzip.open(path_rfam_seed_gz, 'rt', encoding='cp850') as f:
    for line in f:
        if line.strip() and not line.startswith('#') and not line.startswith('//'):
            seq_id = line.split()[0]

            if seq_id not in ids_set:
                ids_set.add(seq_id)
            else:
                if seq_id not in dupId_to_freq_dict:
                    dupId_to_freq_dict[seq_id] = 1
                    dupId_to_suffix_dict[seq_id] = 0
                dupId_to_freq_dict[seq_id] += 1


SS_cons_WUSS_list = []
GC_RF_list = []

# Divide STOCKHOLM files
with gzip.open(path_rfam_seed_gz, 'rt', encoding='cp850') as f:
    for rf_stockholm in [x for x in f.read().split('\n//\n') if x.strip()]:
        to_write = ''
        rf_id = ''
        SS_cons_WUSS = ''
        GC_RF = ''
        for line in rf_stockholm.split('\n'):
            if line.startswith('#=GF AC'):
                # #=GF AC   RF00001
                rf_id = line.split()[2]

            if line.strip() and not line.startswith('#') and not line.startswith('//'):
                seq_id = line.split()[0]

                if seq_id in dupId_to_freq_dict:
                    new_seq_id = seq_id + '_' + str(dupId_to_suffix_dict[seq_id])
                    dupId_to_suffix_dict[seq_id] += 1
                    line = line.replace(seq_id, new_seq_id)
                    print(rf_id, seq_id, new_seq_id)
            elif line.startswith('#=GC SS_cons'):
                # #=GC SS_cons                 :<<<<.<.........
                SS_cons_WUSS = line.strip('\n').split()[2]
            #elif line.startswith('#=GC RF'):
            #    # #=GC SS_cons                 :<<<<.<.........
            #    GC_RF = line.strip('\n').split()[2]

            to_write += line + '\n'

        if rf_id and SS_cons_WUSS:
            # print(rf_id)
            path_rf_stockholm = os.path.join(dir_rfam_stockholm, f'{rf_id}.txt')
            if not os.path.exists(path_rf_stockholm):
                with open(path_rf_stockholm, 'w') as fw:
                    fw.write(to_write + '\n//')

            SS_cons_WUSS_list.append([rf_id, SS_cons_WUSS])

            #if GC_RF:
            #    GC_RF_list.append([rf_id, GC_RF])
        else:
            print('ERRORE', rf_id)

if not os.path.exists(path_SS_cons_WUSS):
    with open(path_SS_cons_WUSS, 'w') as fw:
        fw.write('\n'.join(['\t'.join(x) for x in SS_cons_WUSS_list]) + '\n')

#with open(path_GC_RF, 'w') as fw:
#    fw.write('\n'.join(['\t'.join(x) for x in GC_RF_list]) + '\n')


# Prepare constraints
consensus_symbols = ["A", "U", "C", "G"]
ss_symbols = ["<", ">", "(", ")", "[", "]", "{", "}", "-", "_", ":", ".", ","]
constraint_symbols = ["<", ">", "<", ">", "<", ">", "<", ">", "x", "x", "x", ".", "."]


def cleanConst(const, seq, pairs):
    const = list(const)
    for i, s in enumerate(seq):
        if s == "-" and (const[i] == ">" or const[i] == "<"):
            const[pairs[i]] = "."  # il nucleotide matchato con quello mancante è incognito
            const[i] = "."
    return "".join(const)


def pairedConstraint(const):
    n = 0
    parenthesis = []
    for c in const:
        if c == "<":
            n += 1
            parenthesis.append(n)
        elif c == ">":
            parenthesis.append(n)
            n -= 1
        else:
            parenthesis.append(0)
    pairs = [0] * len(parenthesis)
    for i, pi in enumerate(parenthesis):
        for j, pj in enumerate(parenthesis):
            if pi == pj and i != j and pi != 0:
                pairs[i] = j
                pairs[j] = i
                parenthesis[i] = 0
                parenthesis[j] = 0
                break
    return pairs


def applyConstraint(constraint, sequences):
    finalSequences = {}
    pairs = pairedConstraint(constraint)
    for seq in sequences:
        finalSeq = ""
        finalConst = ""
        if sequences[seq] != "":
            for s, c in zip(sequences[seq], cleanConst(constraint, sequences[seq], pairs)):
                if s != "-":
                    finalSeq += s
                    finalConst += c
            finalSequences[seq] = [finalSeq, finalConst]
    return finalSequences


def convertConstraint(constraint):
    return "".join([constraint_symbols[ss_symbols.index(symbol)] for symbol in constraint])


def createConstraint(ss, consensus):
    constraint = ""
    for s, c in zip(ss, consensus):
        if c in consensus_symbols:
            if s in ss_symbols:
                constraint += s
            else:
                constraint += "-"
        else:
            constraint += "."
    return convertConstraint(constraint)


def constraintFamily(path_rfam_stockholm):
    file = open(path_rfam_stockholm, encoding='cp850').readlines()
    # print(file)

    sequences = {}
    consensus = ''
    ss = ''
    with open(path_rfam_stockholm, encoding='cp850') as f:
        for line in f:
            line = line.strip()
            if line:
                if not line.startswith('#') and not line.startswith('//'):
                    seq_id, seq = line.split()  # To trigger exception in case of problems
                    sequences[seq_id] = seq
                elif line.startswith('#=GC SS_cons'):
                    ss = line.strip('\n').split()[2]
                elif line.startswith('#=GC RF'):
                    consensus = line.strip('\n').split()[2]
    if consensus:
        constraint = createConstraint(ss, consensus)
    else:
        constraint = '.' * len(ss)
    return applyConstraint(constraint, sequences)


def writeConstraint(familyId, dir_constraint, finalSequences):
    path_constraint = os.path.join(dir_constraint, familyId + ".constraint.txt")
    if not os.path.exists(path_constraint):
        with open(path_constraint, "w") as outputFile:
            for seq in finalSequences:
                if seq != "":
                    outputFile.write(">" + seq + "\n" + finalSequences[seq][0] + "\n" + finalSequences[seq][1] + "\n")


for filename in sorted(os.listdir(dir_rfam_stockholm)):
    # print(filename)
    rf_id = filename.split('.')[0]
    writeConstraint(rf_id, dir_rfam_constraint, constraintFamily(os.path.join(dir_rfam_stockholm, filename)))


# RNAfold
rfam_constraint_list = sorted(os.listdir(dir_rfam_constraint))

num_families_div_20 = len(rfam_constraint_list) // 20

print('Folding...')

for num_fam, filename in enumerate(rfam_constraint_list):
    if num_fam % num_families_div_20 == 0:
        print('{:.2f}% ({}/{})'.format((num_fam / len(rfam_constraint_list)) * 100.0, num_fam, len(rfam_constraint_list)))

    rf_id = filename.split('.')[0]

    path_rfam_bear = os.path.join(dir_rfam_bear, rf_id + '.folded.fastb')
    if not os.path.exists(path_rfam_bear):
        path_rfam_folded = os.path.join(dir_rfam_folded, rf_id + '.folded.txt')
        if not os.path.exists(path_rfam_folded):
            path_rfam_constraint = os.path.join(dir_rfam_constraint, filename)

            # Calculate dotbracket
            rna_fold_output = subprocess.check_output(
                [path_RNAfold, '--noPS', path_rfam_constraint],
                ###[path_RNAfold, '--noPS', '-C', path_rfam_constraint],
                # To get strings
                universal_newlines=True
            )

            rf_seq_db = ''
            # Remove the energies
            for x in rna_fold_output.strip('\n>').split('\n>'):
                x_list = x.split('\n')
                x_list[2] = re.sub(r' \((.*?)\)', '', x_list[2])
                rf_seq_db += '>' + '\n'.join(x_list) + '\n'

            # Save folded sequences
            with open(path_rfam_folded, 'w') as fw:
                fw.write(rf_seq_db)

        # Calculate BEAR
        subprocess.call(
            ['java', '-jar', path_BEAR, path_rfam_folded, path_rfam_bear]
        )


# Sequences (without gaps)
for filename in sorted(os.listdir(dir_rfam_bear)):
    # print(filename)
    rf_id = filename.split('.')[0]

    with open(os.path.join(dir_rfam_sequences, rf_id + '.fasta'), 'w') as fw:
        with open(os.path.join(dir_rfam_bear, filename)) as f:
            for line in f:
                if line.startswith('>'):
                    fw.write(line)
                    line = f.readline()
                    fw.write(line)


# Gapped info
rfamId_to_seqId_to_info_to_str_dict = {}

for filename in sorted(os.listdir(dir_rfam_stockholm)):
    # print(filename)
    rf_id = filename.split('.')[0]

    seqId_to_seqGapped = {}

    with open(os.path.join(dir_rfam_stockholm, filename)) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('//'):
                seq_id, seq_gapped = line.split()

                seqId_to_seqGapped[seq_id] = seq_gapped.strip()

    rfamId_to_seqId_to_info_to_str_dict[rf_id] = {}

    path_rfam_bear = os.path.join(dir_rfam_bear, rf_id + '.folded.fastb')
    if not os.path.exists(path_rfam_bear):
        continue

    with open(path_rfam_bear) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip('\n')[1:]
                sequence = f.readline().strip('\n')
                db = f.readline().strip('\n')
                bear = f.readline().strip('\n')

                rfamId_to_seqId_to_info_to_str_dict[rf_id]['>' + seq_id] = {
                    'sequence': seqId_to_seqGapped[seq_id],
                    'db': '',
                    'bear': ''
                }
                if seqId_to_seqGapped[seq_id] == function_mbr.distributeGaps(seqId_to_seqGapped[seq_id], sequence):
                    rfamId_to_seqId_to_info_to_str_dict[rf_id]['>' + seq_id]['db'] += function_mbr.distributeGaps(
                        seqId_to_seqGapped[seq_id], db
                    )
                    rfamId_to_seqId_to_info_to_str_dict[rf_id]['>' + seq_id]['bear'] += function_mbr.distributeGaps(
                        seqId_to_seqGapped[seq_id], bear
                    )
                else:
                    print('ERRORE', rf_id)

pickle.dump(rfamId_to_seqId_to_info_to_str_dict, gzip.open(path_gapped_fam_dict_pickle, 'wb'))
