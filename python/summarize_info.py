file1 = 'sra_result.csv'
file2 = 'SraRunInfo.csv'

tmp2 = {}
with open(file2, 'r') as f:
	next(f)
	for line in f:
		line = line.strip('\n').split(',')
		run = line[0]
		sample = line[24]
		tmp2[sample] = run

organs = ['Adipose', 'Heart', 'Kidney', 'Liver']
tmp1 = {k:[] for k in organs}
with open(file1, 'r') as f:
	next(f)
	for line in f: 
		line = line.strip('\n').split(',')
		experiment_title = line[1].replace('"','').replace(' ','')
		sample = line[7].replace('"','')
		run = tmp2[sample]
		
		exp_organ = None
		for o in organs:
			if o in experiment_title:
				exp_organ = o
				break

		tmp1[exp_organ].append((experiment_title, run))

with open('summary_info.txt', 'w') as f:
	for o, value in tmp1.items():
		for c, (run, info) in enumerate(value[::-1]):
			f.write(f'{o}-{c+1},{info},{run}\n')