#!/usr/bin/env python

from Bio import SeqIO
import argparse

def parse_cli():
    parser = argparse.ArgumentParser(description="Concatenate fasta formatted alignments")

    parser.add_argument('-i', type=str, dest='files', nargs='+',
                        help='names of 2 or more files to be concatenated')

    parser.add_argument('-o', type=str, dest='outfile',
                        help='output file')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    opts = parse_cli()

    superaln = {}

# NOTE: it looks like there are some functions in Bio.Align
# that can accomplish this task...

    for file in opts.files:
        for seq in SeqIO.parse(file, 'fasta'):
            if seq.id not in superaln:

                superaln[seq.id] = seq
            else:
                superaln[seq.id] += seq

    SeqIO.write(superaln.values(), opts.outfile, 'fasta')