import os
import csv
import pdb
from Bio import SeqIO
from sys import argv
import argparse


def parse_cli():
    parser = argparse.ArgumentParser(description="Concatenate fasta formatted alignments")

    parser.add_argument('-i', type=str, dest='file',
                        help='fasta file to be renamed')

    parser.add_argument('-s', type=str, dest='species_list',
                        help='output file')

    parser.add_argument('-o', type=str, dest='outfile',
                        help='output file')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    opts = parse_cli()

    names = {}
    out_seqs = []

    for row in csv.reader(open(opts.species_list), delimiter=','):
        if 'letter code' in  row[0]:
            continue

        sp = row[1].replace(' ', '_')
        names[row[0]] = sp + '_' +  row[2]

    for seq in SeqIO.parse(opts.file, 'fasta'):
            new_name = names[seq.id]
            seq.description = seq.id
            seq.id = new_name
            out_seqs.append(seq)

    SeqIO.write(out_seqs, opts.outfile, 'fasta')