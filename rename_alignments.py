import os
import csv
import pdb
from Bio import AlignIO
from sys import argv
from Bio.Alphabet import ProteinAlphabet
in_file = 'species_list.csv'

names = {}
for row in csv.reader(open(in_file), delimiter=','):
    if 'letter code' in  row[0]:
        continue

    sp = row[1].replace(' ', '_')
    cls = row[2]
    names[row[0]] = sp

work_dir = 'align-SH'


for file in os.listdir(work_dir):
    if file.endswith('aligned.fasta'):
        #outseqs = []
        alpha = ProteinAlphabet()
        alignment = AlignIO.read(work_dir+'/'+file, 'fasta', alphabet=alpha)
        for seq in alignment:
            seq.id = names[seq.id]
            #outseqs.append(seq)
            #pdb.set_trace()
        outfile = work_dir + '/' + file + '.renamed.nex'
        #pdb.set_trace()
        outfile_h = open(outfile, 'w')
        AlignIO.write([alignment], outfile_h, 'nexus')
        outfile_h.close()