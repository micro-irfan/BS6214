from collections import defaultdict

file = 'LiverTop100PCA.txt'
liverTop = set()
liverRank = []
with open('LiverTop100PCA.GO.txt', 'w') as write_file:
    with open(file, 'r') as f:
        for line in f:
            line = line.strip('\n')
            liverTop.add(line)
            liverRank.append(line)
            write_file.write(line.split('.')[0] + '\n')

file = 'pharmacogenes.long.txt'
pharmacogenes = set()
with open(file, 'r') as f:
    for line in f:
        line = line.strip('\n')
        pharmacogenes.add(line)

file = 'allGeneLength.txt'
geneDict = defaultdict(list)
count = 0
with open(file, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        gene = line[0]
        if gene in liverTop: 
            # if line[1] in geneDict.keys():
            #     # print ('REPEAT', gene, line[0], geneDict[gene])
            #     if gene == 'GPI': continue

            geneDict[line[1]].append(gene)
            if line[1] in pharmacogenes:
                index = liverRank.index(line[0])
                print ('Pharmacogenes!', index, line[1], line[0], geneDict[line[1]])
                count += 1
            
print (count)
for g, v in geneDict.items():
    print (g, v)

