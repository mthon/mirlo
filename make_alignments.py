#!/usr/bin/python
import os
import pdb
import csv
import argparse

from Bio import SeqIO

# seq_file = 'all.fasta'
# clust_file = 'all.clusters'
# in_file = 'species_list.csv'
# out_dir = 'alignments'


def parse_cli():
    parser = argparse.ArgumentParser(description="Rename sequences")

    parser.add_argument('-i', type=str, dest='seq_file', default='all.fasta',
                        help='protein sequences in fasta format')

    parser.add_argument('-c', type=str, dest='clust_file', default='all.clusters',
                        help='protein clusters file from MCL')

    parser.add_argument('-s', type=str, dest='species_list',
                        default='species_list.csv',
                        help='text file containing the file names and three letter codes')

    parser.add_argument('-o', type=str, dest='out_dir',
                        default='alignments',
                        help='output directory')

    opts = parser.parse_args()

    return opts


def do_it(opts):

    try:
        os.stat(opts.out_dir)
    except:
        os.mkdir(opts.out_dir)

    num_spp = 0
    for row in csv.reader(open(opts.species_list), delimiter=','):
        if 'letter code' in row[0]:
            continue
        if len(row[0]) > 0:
            num_spp += 1

    file_list = []
    seqs = {}
    #pdb.set_trace()

    for seq in SeqIO.parse(open(opts.seq_file), 'fasta'):
        seqs[seq.id] = seq

    clust_id = 0

    for row in csv.reader(open(opts.clust_file), delimiter='\t'):
        clust_id += 1
        counts = {}
        include = True
        for seq in row:
            spid = seq[0:seq.index('::')]
            if spid not in counts:
                counts[spid] = 0
            counts[spid] = counts[spid] + 1

        for spid, count in counts.iteritems():
            if count > 1:
                include = False
        if len(counts.keys()) < num_spp:
            include = False

        if include:
            print 'saving cluster %s' % str(clust_id)
            out_file = '%s/%s.fasta' % (opts.out_dir, str(clust_id))
            out_seqs = []

            for seqid in sorted(row):
                seq = seqs[seqid]
                seq.id = seqid[0:seqid.index('::')]
                out_seqs.append(seq)
            SeqIO.write(out_seqs, out_file, 'fasta')
            file_list.append(out_file)

    print 'starting alignments'

    #pdb.set_trace()
    for file in file_list:

        print 'aligning file %s' % file
        out_file = file + '.aligned.fasta'
        command = 'mafft --amino --auto  --inputorder %s > %s' % (file, out_file)
        os.system(command)

if __name__ == '__main__':
    opts = parse_cli()
    do_it(opts)
