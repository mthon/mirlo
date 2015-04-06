import os
import csv
import pdb
from Bio import AlignIO
from sys import argv
from Bio.Alphabet import ProteinAlphabet

import argparse

def parse_cli():
    parser = argparse.ArgumentParser(description="Rename alignments")

    parser.add_argument('-i', type=str, dest='work_dir', default='alignments',
                        help='directory containing aligned sequences')

    parser.add_argument('-s', type=str, dest='in_file',
                        help='text file containing the file names and three letter codes')

    opts = parser.parse_args()

    return opts

if __name__ == '__main__':
    cli_args = parse_cli()
    names = {}
    for row in csv.reader(open(cli_args.in_file), delimiter=','):
        if 'letter code' in  row[0]:
            continue

        sp = row[1].replace(' ', '_')
        cls = row[2]
        names[row[0]] = sp

    for file in os.listdir(cli_args.work_dir):
        if file.endswith('aligned.fasta'):
            #outseqs = []
            alpha = ProteinAlphabet()
            alignment = AlignIO.read(cli_args.work_dir+'/'+file, 'fasta', alphabet=alpha)
            for seq in alignment:
                seq.id = names[seq.id]
                #outseqs.append(seq)
                #pdb.set_trace()
            outfile = cli_args.work_dir + '/' + file + '.renamed.nex'
            #pdb.set_trace()
            outfile_h = open(outfile, 'w')
            AlignIO.write([alignment], outfile_h, 'nexus')
            outfile_h.close()