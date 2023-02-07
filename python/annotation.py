
header = 'SRR2087282.out.bam	SRR2087283.out.bam	SRR2087284.out.bam	SRR2087285.out.bam	SRR2087286.out.bam	SRR2087287.out.bam	SRR2087288.out.bam	SRR2087289.out.bam	SRR2087290.out.bam	SRR2087291.out.bam	SRR2087292.out.bam	SRR2087293.out.bam	SRR2087294.out.bam	SRR2087295.out.bam	SRR2087296.out.bam	SRR2087297.out.bam	SRR2087298.out.bam	SRR2087299.out.bam	SRR2087300.out.bam	SRR2087301.out.bam	SRR2087302.out.bam	SRR2087303.out.bam	SRR2087304.out.bam	SRR2087305.out.bam	SRR2087306.out.bam	SRR2087307.out.bam	SRR2087308.out.bam	SRR2087309.out.bam	SRR2087310.out.bam	SRR2087311.out.bam	SRR2087312.out.bam	SRR2087313.out.bam	SRR2087314.out.bam	SRR2087315.out.bam	SRR2087316.out.bam	SRR2087317.out.bam	SRR2087318.out.bam	SRR2087319.out.bam	SRR2087320.out.bam	SRR2087321.out.bam	SRR2087322.out.bam	SRR2087323.out.bam	SRR2087324.out.bam	SRR2087325.out.bam	SRR2087327.out.bam	SRR2087328.out.bam	SRR2087329.out.bam	SRR2087330.out.bam	SRR2087331.out.bam	SRR2087332.out.bam	SRR2087333.out.bam	SRR2087334.out.bam	SRR2087335.out.bam	SRR2087336.out.bam	SRR2087337.out.bam	SRR2087338.out.bam	SRR2087339.out.bam	SRR2087340.out.bam	SRR2087341.out.bam	SRR2087342.out.bam	SRR2087343.out.bam	SRR2087344.out.bam	SRR2087345.out.bam	SRR2087346.out.bam	SRR2087347.out.bam	SRR2087348.out.bam	SRR2087349.out.bam	SRR2087350.out.bam	SRR2087351.out.bam	SRR2087352.out.bam	SRR2087353.out.bam	SRR2087354.out.bam	SRR2087355.out.bam	SRR2087356.out.bam	SRR2087357.out.bam	SRR2087358.out.bam	SRR2087359.out.bam	SRR2087360.out.bam	SRR2087361.out.bam	SRR2087362.out.bam	SRR2087363.out.bam	SRR2087364.out.bam	SRR2087365.out.bam	SRR2087366.out.bam	SRR2087367.out.bam	SRR2087368.out.bam	SRR2087369.out.bam	SRR2087370.out.bam	SRR2087371.out.bam	SRR2087372.out.bam	SRR2087373.out.bam	SRR2087374.out.bam	SRR2087375.out.bam'
header = header.split('\t')


file = 'summary_info.txt'
with open(file, 'r') as f, open('annotationR.csv', 'w') as write_file:
	tmp = {}
	for line in f:
		line = line.strip('\n').split(',')
		expriment = line[0].split('-')[0]
		accession = line[1]
		tmp[accession] = expriment

	write_file.write(f'run,expriment\n')
	for h in header:
		accession = h.replace('.out.bam','')
		expriment = tmp[accession]
		write_file.write(f'{h},{expriment}\n')

