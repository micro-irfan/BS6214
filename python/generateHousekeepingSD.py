import math
import pandas as pd 
import numpy as np
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt

tpm = pd.read_csv('Human_genes.tpm.csv') 
tpm.rename(columns={"Unnamed: 0": "geneID"}, inplace=True)
geneName = pd.read_csv('allGeneLength.txt', sep='\t') 
geneName.drop(columns=['start','end','length'], inplace=True)
tpm = pd.merge(tpm, geneName, how='inner', on = 'geneID')
tpm.set_index(['geneName', 'geneID'], inplace=True)

## Remove rows sum == 0
tpm = tpm.loc[~(tpm==0).all(axis=1)]
tpm = tpm.T

file = 'annotationR.csv'
rename_dict = {}
with open(file, 'r') as f:
    next(f)
    for line in f:
        line = line.strip('\n').split(',')
        rename_dict[line[0]] = f'{line[1]}-{line[2]}'

tpm1 = tpm.T.rename(columns=rename_dict)

housekeeping_genes1 = ['EMC7', 'C1orf43', 'CHMP2A', 'GPI', 'PSMB2', 'PSMB4', 'RAB7A', 'REEP5', 'SNRPD3', 'VCP', 'VPS29']
organs = ['Liver', 'Kidney', 'Adipose', 'Heart']

for c_genes, housekeeping_genes in enumerate([housekeeping_genes1]):
    housekeeping_values = {}
    for g in housekeeping_genes:
        try:
            i = tpm1.index.get_loc(g)
            housekeeping_values[g] = tpm1.iloc[i]
        except:
            print (f'{g} Not found')
    
    housekeeping_genes_values = {}
    for g in housekeeping_genes:
        tmp_df = housekeeping_values[g]
        tmp = defaultdict(list)
        for i in tmp_df.columns:
            value = np.log2(tmp_df[i][0])
            organ = i.split('-')[0]
            tmp[organ].append(value)
        housekeeping_genes_values[g] = tmp

    sd_dict = {}
    for g, values in housekeeping_genes_values.items():
        tmp = {}
        for organ, v in values.items():
            median = np.median(v)
            sd = np.std(v)
            tmp[organ] = sd
        sd_dict[g] = tmp

    data = [[0 for i in range(len(housekeeping_genes))] for i in range(4)]

    for c1, (k,v) in enumerate(sd_dict.items()):
        for o, v1 in v.items():
            c2 = organs.index(o)
            data[c2][c1] = round(v1, 5)
            #print (f'{k}:{o}:{round(v1, 3)}')

    df = pd.DataFrame(data, columns = housekeeping_genes)
    df['Organ'] = organs
    df.set_index('Organ', inplace = True)
    df.to_csv(f'Housekeeping-genes.{c_genes}.csv', index=True)