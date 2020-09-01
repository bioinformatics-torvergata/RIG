from collections import Counter
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import subprocess

#funtion to run BlustClust (filter on identity score)
#folder = seq_str_families/sequence/
#identity = threshold expressed by % (ex. 50)
#return a folder with selected sequences
def run_blustClust(folder, identity, name):
    if 'not_similar_'+name+'_'+identity not in os.listdir('./'):
        os.mkdir('not_similar_'+name+'_'+identity+'/')
    os.chdir(folder)
    os.system('for file in *; do blastclust -i $file -o ../not_similar_'+name+'_'+identity+'/$file -p F -S '+identity+' ; done')
    os.chdir('../')
    print('BlustClust DONE!')
    

#function to filter for number of sequences in an RNA family after BlustClust
#folder = folder generated from run_blustClust
#n_seq_family = minimum number of sequences per family
#identity = threshold expressed by % (ex. 50)
def filter_n_seq(folder, n_seq_family, name, identity):
    if 'filter_n_seq_'+name+'_'+identity not in os.listdir('./'):
        os.mkdir('filter_n_seq_'+name+'_'+identity+'/')
    list_fam=os.listdir(folder)
    #New folder with only families with more than n_seq_family members
    for fam in list_fam:
        output=subprocess.check_output("wc -l "+folder+fam, shell=True)
        if int(output.strip().split()[0])>=int(n_seq_family):
            os.system('cp '+folder+'/'+fam+' filter_n_seq_'+name+'_'+identity+'/'+fam)
    print ('Filter Nseq DONE!')
    
    
  
#function to select the bear sequences after BlustCLust filtering
def get_bear(folder, folder_bear, name, identity):    
    #folder = folder returned by filter_n_seq
    #folder_bear = folder seq_str_families/bear
    if 'bear_filtered_'+name+'_'+identity not in os.listdir('./'):
        os.mkdir('bear_filtered_'+name+'_'+identity+'/')
    list_fam_filter=os.listdir(folder)
    for fam_clean in list_fam_filter:
        seq=[]
        f=open(folder+fam_clean).readlines()
        for line in f:
            seq.append(line.split()[0])

        o=open('bear_filtered_'+name+'_'+identity+'/'+fam_clean, "w")
        f2=open(folder_bear+fam_clean)
        line=f2.readline()
        while(line):
            if line[0]==">" and line[1:-1] in seq:
                o.write(line)
                line=f2.readline()
                o.write(line)
                line=f2.readline()
                o.write(line)
                line=f2.readline()
                o.write(line)
                line=f2.readline()
            else:
                line=f2.readline()
        o.close()
        
    print ('Bear sequence research DONE!')
    

###MARCO####
def distributeGaps(gappedReference, ungappedString):
    assert len(gappedReference.replace('-','')) == len(ungappedString), 'ungapped strings should be equal' 
    result = list(ungappedString)
    gaplist = [ m.start() for m in re.finditer('-', gappedReference)]

    for gap in gaplist:    
        result.insert(gap, '-')
    result = "".join(result)    
    return result

assert(distributeGaps('--abcdef-g-', 'bombasi') == '--bombas-i-')


#Function to Add gap based on RFAM alignment from seed
#folder=bear_filtered
#seed_rfam
def add_gap(folder, seed_rfam, name, identity):
    list_fam_filter2=os.listdir(folder)
    if 'bear_alignment_'+name+'_'+identity not in os.listdir('./'):
        os.mkdir('bear_alignment_'+name+'_'+identity+'/')
    c=0
    for fam in list_fam_filter2:
        c+=1
        o=open('bear_alignment_'+name+'_'+identity+'/'+fam, "w")
        f=open(folder+fam)
        line=f.readline()
        while(line):
            if line[0]==">":
                seq_name=line[1:-1]
                #output=subprocess.check_output("wc -l not_similar/"+fam, shell=True)
                seq_alignment=subprocess.check_output('grep '+seq_name+' '+seed_rfam, shell=True)
                line=f.readline()
                test=distributeGaps(seq_alignment.split()[1], line[0:-1])
                line=f.readline()
                line=f.readline()
                bear_seq=line
                gap_pos=[]
                for i, el in enumerate(test):
                    if el=='-':
                        gap_pos.append(i)
                for i in gap_pos:
                    bear_seq=bear_seq[:i] + "-" + bear_seq[i:]
                o.write(bear_seq)
            else:
                line=f.readline()
        o.close()

        #print fam+"\t"+str(c)+" famiglie su "+str(len(lista_fam_filter2))
    
    print ('Alignment with gap DONE!')
    print ('Bear sequences with gap in folder bear_alignment/')
    



###MARCO###
def decode(bear):
    alph_bear={'abc':'a', 'def':'A', 'ghi=':'=',
          'lmnop':'l', 'qrstu':'L', 'vwxyz^':'^',
          '!"#': 'i', '$%&':'I', '\'()+':'+',
          '234':'n', '567':'N', '890>':'>',
          'ABC':'s', 'DEF':'S', 'GHIJ':'~',
          'KLMN':'b', 'OPQR':'B', 'STUVW':'|',
          'YZ~':'y', '?_|':'Y', '/\\@':'@',
          '{}[]':'[', ':':':' , '-':'-'}
    
    result=""
    for ch in bear:
        for key in alph_bear:
            if ch in key:
                result += alph_bear[key]
    return result


#deconding2 from text (same as decode but from text file)
def decode_from_file(bear, file_name):
    f=open(file_name).readlines()
    alph_bear=[i.split() for i in f]
    alph_bear.append(['-','-'])
    result=""
    for ch in bear:
        for group in alph_bear:
            if ch in group[0]:
                result += group[1]
    return result


#convert from bear to New Alphabet
#folder = bear_alignment(bear_alignment/bear_alignment___62)
def convert_new_bear(folder, name, identity):
    fam_bear=os.listdir(folder)
    if 'bear_new_alignment_'+name+'_'+identity not in os.listdir('./'):
        os.mkdir('bear_new_alignment_'+name+'_'+identity+'/')
    for fam in fam_bear:
        o=open('bear_new_alignment_'+name+'_'+identity+'/'+fam, "w")
        f=open(folder+fam)
        line=f.readline()
        while(line):
            o.write(decode(line)+"\n")
            line=f.readline()
        o.close()
    print ('Decoding BEAR DONE!')

#convert from bear to New Alphabet (FILE)
#folder = bear_alignment(bear_alignment/bear_alignment___62)
#file_name = alph_mapping.tsv
def convert_new_bear_file(folder, file_name, name, identity):
    fam_bear=os.listdir(folder)
    if 'bear_new_alignment_'+name+'_'+identity not in os.listdir('./'):
        os.mkdir('bear_new_alignment_'+name+'_'+identity+'/')
    for fam in fam_bear:
        o=open('bear_new_alignment_'+name+'_'+identity+'/'+fam, "w")
        f=open(folder+fam)
        line=f.readline()
        while(line):
            o.write(decode_from_file(line, file_name)+"\n")
            line=f.readline()
        o.close()
    print ('Decoding BEAR from File DONE!')
    

#Create Blocks from bear alignment
#folder = bear_alignment/bear_new_alignment_$alph_$id/
def make_blocks(folder, name, identity):
    if 'blocks_new_bear_'+name+'_'+identity not in os.listdir('./'):
        os.mkdir('blocks_new_bear_'+name+'_'+identity+'/')
    families=os.listdir(folder)
    for fam in families:
        o=open('blocks_new_bear_'+name+'_'+identity+'/'+fam, "w")
        f=open(folder+fam).readlines()
        zipped=zip(*f)
        v=[]
        for col in zipped:
            if '-' not in col:
                c=Counter(col)
                for key in c:
                    if float(c[key])/float(len(col))>0.5:
                        v.append(col)
        new=zip(*v)
        for seq in new:
            o.write(''.join(seq))
            

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
                
    fr_expected.to_csv('expected_frequencies.tsv',sep="\t")
    return fr_expected
    print ('Expected_frequencies DONE!')
    
#observed_substitution(f_ij)
#folder=Blocks/blocks_new_bear_$alph_$id
#v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']

def observed_substitution(folder, v_bear, name, identity):
    substitution=pd.DataFrame(1.0, columns=v_bear, index=v_bear)
    list_=os.listdir(folder)
    for fam in list_:
        #print fam
        f=open(folder+fam).readlines()
        v1, v2 = np.triu_indices(len(f))
        for k in range(0, len(v1)):
            if v1[k]!=v2[k]:
                for i, el in enumerate(f[v1[k]].strip()):
                    substitution.ix[el, f[v2[k]].strip()[i]]+=1.0
                    substitution.ix[f[v2[k]].strip()[i], el]+=1.0
    substitution.to_csv('substitution_'+name+'_'+identity+'.tsv', sep="\t") 
    return substitution
    print ('Substitution matrix DONE!')

 
 
 #substitution=dataframe returned by observed_substitution
 #v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
def make_q(substitution, v_bear, name, identity):
    number_couple=0
    for j in range(0, len(v_bear)):
        number_couple+=substitution.iloc[j,j:].sum()
    #print number_couple
    q_ij=substitution.divide(number_couple)
    q_ij.to_csv('q_ij_'+name+'_'+identity+'.tsv', sep="\t")
    return q_ij
    print ('q_ij DONE!')
    

#q_ij returned by make_q
#v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
def make_p(q_ij, v_bear):
    p_i=dict.fromkeys(v_bear, 0)
    for char in v_bear:
        same=q_ij[char][char]
        others=(sum(q_ij[char])-same)/2
        p_i[char]+=same+others
    return p_i
    print ('p_i dict DONE!')

#p_i returned by make_p
#v_bear=['a','A','=','l','L','^','i','I','+','n','N','>','s','S','~','b','B','|','y','Y','@','[', ':']
def make_e(p_i, v_bear, name, identity):
    e_ij=pd.DataFrame(1.0, columns=v_bear, index=v_bear)
    for char in v_bear:
        for char2 in v_bear:
            if char==char2:
                e_ij.ix[char, char2]=p_i[char]*p_i[char2]
            else:
                e_ij.ix[char, char2]=2*p_i[char]*p_i[char2]
    
    e_ij.to_csv('E_ij_'+name+'_'+identity+'.tsv', sep='\t')
    return e_ij
    
    print ('e_ij DONE!')
    
    
#Make Matrix
#freq_observed=dataframe of observed frequencies (q_ij)
#fr_expected=dataframe of expected frequencies (e_ij)
def make_matrix(freq_observed, fr_expect, name, identity):
    #Probability maytrix
    probability_matrix=freq_observed.divide(fr_expect)
    probability_matrix.to_csv('probability_matrix_'+name+'_'+identity+'.tsv', sep="\t")
    print ('Probability matrix DONE!')
    
    #Score Matrix
    mbr_new=probability_matrix.applymap(np.log2)
    mbr_new.to_csv('MBR_'+name+'_'+identity+'.tsv', sep="\t")
    print ('MBR DONE!')
    return mbr_new
    
    
def Expected_score(s_ij, p_i):
    E=0
    for i, char in enumerate(s_ij.columns.values):
        for i2, char2 in enumerate(s_ij.columns.values):
            if i2>=i:
                #print char, char2
                E+=s_ij.ix[char, char2]*p_i[char]*p_i[char2]
    #print E
    return E


def entropy(q_ij, s_ij):
    H=0
    length=len(q_ij.columns.value_counts())
    for j in range(0, length):
        H+=(q_ij.iloc[j,j:]*s_ij.iloc[j,j:]).sum()
    #print H
    return H
    
    
def make_heatmap(S_ij, name, identity):
    sns.set(font_scale=2.0)
    plt.figure(figsize=(10,10))
    ax=sns.heatmap(S_ij, xticklabels=1, yticklabels=1,cmap="YlGnBu")
    plt.savefig('matrix_'+name+'_'+identity+'.pdf')
    #plt.show()
    plt.close()
    

#Function that start from sequence/structure of all RFAM families and create a new MBR from alphabet file
#folder = folder with RFAM RNA sequences
#folder_bear = folder with RFAM RNA structure in bear
#RFAM_seed_file = Rfam seed downloaded from RFAM database
#id_blustClust = Sequence identity for filter
#filter_nSeq = threshold on number of sequences in a family after filtering
#file_alph = alphabet file with bear mapping
#file_info = output file with information about the MBR

def BlustClust_filter_alignment(folder, folder_bear, RFAM_seed_file, id_blustClust, filter_nSeq, file_alph):
    name=file_alph.split('.')[0]
    run_blustClust(folder, id_blustClust, name)
    filter_n_seq('not_similar_'+name+'_'+id_blustClust+'/', filter_nSeq, name, id_blustClust)
    get_bear('filter_n_seq_'+name+'_'+id_blustClust+'/', folder_bear, name, id_blustClust)
    add_gap('bear_filtered_'+name+'_'+id_blustClust+'/', RFAM_seed_file, name, id_blustClust)
    convert_new_bear_file('bear_alignment___'+id_blustClust+'/', file_alph, name, id_blustClust)


def Make_MBR_from_blocks(blocks_folder, id_blustClust, file_alph, file_info):
    name=file_alph.split('.')[0]
    f=open(file_alph).readlines()
    o=open(file_info, "w")
    #blocks_folder='blocks_new_bear_'+name+'_'+id_blustClust+'/'
    
    v_char=[x.split()[1] for x in f]
    sost=observed_substitution(blocks_folder, v_char, name, id_blustClust)
    q_ij=make_q(sost, v_char, name, id_blustClust) 
    p_i=make_p(q_ij, v_char)
    e_ij=make_e(p_i, v_char, name, id_blustClust)
    S_ij=make_matrix(q_ij, e_ij, name, id_blustClust)
    E=Expected_score(S_ij, p_i)
    H=entropy(q_ij, S_ij)
    o.write('Expected_score:\t'+str(E)+'\nEntropy:\t'+str(H))
    o.close()
    make_heatmap(S_ij, name, id_blustClust)
    
    