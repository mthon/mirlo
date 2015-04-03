#!/usr/bin/env python
#import os
import sys
import csv
import pdb
import argparse
from Bio import SeqIO

def parse_cli():
    parser = argparse.ArgumentParser(description="Rename sequences")
    parser.add_argument('-d', type=str, dest='in_dir', default='fasta',
        help='directory of fasta files containing the sequences to be renamed')
    parser.add_argument('-s', type=str, dest='species_list',
        default='species_list.csv',
        help='text file containing the file names and three letter codes')
    parser.add_argument('-o', type=str, dest='outfile', default='all.fasta',
        help='output file')

    opts = parser.parse_args()

    return opts

if __name__ == '__main__':

    cli_args = parse_cli()

    all_seqs = []
    for row in csv.reader(open(cli_args.species_list), delimiter=','):
        if 'letter code' in  row[0]:
            continue
        file_h = open(cli_args.in_dir + '/' + row[3])
        for seq in SeqIO.parse(file_h, 'fasta'):
            seq.id = row[0] + '::' + seq.id
            all_seqs.append(seq)

    out_h = open(cli_args.outfile, 'w')
    SeqIO.write(all_seqs, out_h, 'fasta')
