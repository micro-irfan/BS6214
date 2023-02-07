from collections import defaultdict
import math

housekeeping_genes = ['EMC7', 'C1orf43', 'CHMP2A', 'GPI', 'PSMB2', 'PSMB4', 'RAB7A', 'REEP5', 'SNRPD3', 'VCP', 'VPS29']

file = 'allGeneLength.txt'
geneDict = {}
with open(file, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        gene = line[1]
        if gene in housekeeping_genes: 
            if gene in geneDict.keys():
                print (gene, line[0], geneDict[gene])
                if gene == 'GPI': continue
            geneDict[gene] = line[0]

geneDict = {v:k for k,v in geneDict.items()}
print (geneDict)

file = 'annotationR.csv'
annotateDict = defaultdict(list)
with open(file, 'r') as f:
    next(f)
    for line in f: 
        line = line.strip('\n').split(',')
        annotateDict[line[1]].append(line[0])

tissues = list(annotateDict.keys())
matrix = 'Human_genes.mx'
with open(matrix, 'r') as f:
    print_matrix = [[0 for i in range(len(housekeeping_genes))] for j in range(4)]
    for c, line in enumerate(f): 
        line = line.strip('\n').split('\t')
        if c == 0:
            headers = []
            for h in range(1, len(line)):
                for k,v in annotateDict.items():
                    if line[h] in v:
                        headers.append(k)
                        break
            continue

        if not line[0] in geneDict.keys():
            continue

        geneID = line[0]
        values = defaultdict(list)
        for i in range(1, len(line)):
            values[headers[i-1]].append(math.log10(float(line[i])))
        
        whole_geometric_mean = 1
        geometric_mean_dict = {}
        for k, v in values.items():
            multiply = 1
            for i in v:
                multiply *= i 
                whole_geometric_mean *= i

            geometric_mean = (multiply)**(1/len(v))
            geometric_mean_dict[k] = geometric_mean

        whole_geometric_mean = whole_geometric_mean**(1/(len(line)-1))

        print (geometric_mean_dict)
        print (whole_geometric_mean)

        tissues = ['Kidney', 'Heart', 'Adipose', 'Liver']
        
        index = housekeeping_genes.index(geneDict[geneID])
        for k,v in geometric_mean_dict.items():
            tissue_index = tissues.index(k)
            print_matrix[tissue_index][index] = round(v/whole_geometric_mean, 5)

for c, m in enumerate(print_matrix):
    print (tissues[c], [float(f'{i:.3f}') for i in m])