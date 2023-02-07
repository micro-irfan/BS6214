from collections import defaultdict

file = 'pharmacogenes.long.txt'

genes = set()
with open(file, 'r') as f:
    for line in f:
        genes.add(line.strip('\n'))

file = 'allGeneLength.txt'
geneDict = {g:[] for g in genes}
with open(file, 'r') as f:
    next(f)
    for line in f:
        line = line.strip('\n').split('\t')
        gene = line[1]
        if gene in genes: 
            geneDict[gene].append(line[0])
            continue

geneID = set()
with open('pharmacogenes.annotate.txt', 'w') as f:
    for g, value in geneDict.items():
        if not value: 
            f.write(f'{g},NOT FOUND\n')
            print (f'{g} not found!')
            continue

        f.write(f'{g},{";".join(value)}\n')
        for v in value:
            geneID.add(v)


file1 = 'Human_genes.mx'
file2 = 'Human_genes.pharmacogenes.mx'

with open(file1, 'r') as f1:
    with open(file2, 'w') as f2:
        for c, line in enumerate(f1):
            if c == 0:
                f2.write(line)
                continue

            col = line.strip('\n').split('\t') 
            if col[0].replace('"','') in geneID:
                f2.write(line)

file1 = 'Human_genes.tpm.csv'
file2 = 'Human_genes.pharmacogenes.tpm.csv'

with open(file1, 'r') as f1:
    with open(file2, 'w') as f2:
        for c, line in enumerate(f1):
            if c == 0:
                f2.write(line)
                continue
                
            col = line.strip('\n').split(',') 
            if col[0].replace('"','') in geneID:
                f2.write(line)