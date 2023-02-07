
file = 'annotationR.csv'
liver = []
with open(file, 'r') as f:
	next(f)
	for c, line in enumerate(f):
		line = line.strip('\n').split(',')
		if line[1] == 'Liver':
			liver.append(line[0])

file = "Human_genes.mx"
liverFile = 'Human_genes.liver.mx'
with open(file, 'r') as f, open(liverFile, 'w') as write_file:
	for c, line in enumerate(f):
		line = line.strip('\n').split('\t')
		if c == 0: 
			headers = line
			indexes = []
			for i in liver: 
				indexes.append(headers.index(i))	

			row_value = "\t".join(liver)
			write_file.write(f'GeneID\t{row_value}\n')
			continue

		row_value = []
		for i in indexes:
			row_value.append(line[i])

		row_value = "\t".join(row_value)
		write_file.write(f'{line[0]}\t{row_value}\n')