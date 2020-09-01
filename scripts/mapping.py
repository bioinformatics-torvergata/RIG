import sys

# How to run 'python mapping file_alphabet file_bear > out_file'

file_alph = sys.argv[1]
file_bear = sys.argv[2]

# Deconding from text

f = open(file_alph).readlines()
alph_bear = [i.split() for i in f]
alph_bear.append(['-', '-'])

f = open(file_bear)

line = f.readline()
while line:
    print(line.strip())
    line = f.readline()
    print(line.strip())
    line = f.readline()
    print(line.strip())
    line = f.readline()
    bear = line.strip()

    result = ""
    for ch in bear:
        for group in alph_bear:
            if ch in group[0]:
                result += group[1]

    print(result)
    line = f.readline()
