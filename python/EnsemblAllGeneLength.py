#!/usr/bin/env python3
import re
import sys

file = '../BamFiles/genes.gtf'

with open(file, mode='r') as f, open('allGeneLength.txt', mode='w') as write_file:
    write_file.write('geneID\tgeneName\tstart\tend\tlength\n')
    for line in f:
        if not line.startswith('#'): 
            fields = line.split('\t')
            annotations = fields[8].split(' ')
            #print (annotations)
            geneID   = re.sub('[";]', '', annotations[1])
            if fields[2] == 'gene':
                geneName = re.sub('[";]', '', annotations[7])
                write_file.write('%s\t%s\t%s\t%s\t%d\n' %
                      (geneID, geneName, fields[3], fields[4],
                       int( int(fields[4]) - int(fields[3]) )
                      ))