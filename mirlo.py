#!/usr/bin/python
import os
import sys
import csv
import argparse
import subprocess
from Queue import Queue
from threading import Thread
from time import sleep
import StringIO

from Bio import SeqIO
from Bio.AlignIO import convert
from Bio.Phylo import read

import pdb


def parse_cli():
    parser = argparse.ArgumentParser(description="Make alignments of single copy gene families")

    parser.add_argument('-i', type=str, dest='seq_dir', default='fasta',
                        help='directory containing fasta files of proteins')

    parser.add_argument('-c', type=str, dest='clust_file', default='Orthogroups.csv',
                        help='the Orthogroups.csv file from OrthoFinder')

    parser.add_argument('-p', type=str, dest='prottest',
                        help='path to prottest jar file')

    parser.add_argument('-o', type=str, dest='out_dir',
                        default='alignments',
                        help='output directory')

    parser.add_argument('-t', type=str, dest='threads',
                        default='1',
                        help='number of parallel processes')

    opts = parser.parse_args()

    return opts

def parse_fasta(fasta_dir):

    seqs = {}
    num_spp = 0 # test

    if not os.path.isdir(fasta_dir):
        sys.exit('Error: %s is not a directory' % fasta_dir)

    for file in os.listdir(fasta_dir):
        file_path = os.path.join(fasta_dir, file)

        if not os.path.isfile(file_path):
            #sys.exit('Error: cannot read file %s' % file_path)
            continue
        num_spp += 1
        for seq in SeqIO.parse(file_path, 'fasta'):
            seq.source = file
            seq.seq = seq.seq.rstrip('*')
            seqs[seq.id] = seq

    return seqs, num_spp

def parse_clusters(seqs, clust_file, num_spp):
    print 'finding clusters'
    cluster_count = 0
    file_list  = []
    #file_names = []

    for row in csv.reader(open(clust_file), delimiter='\t'):

        if row[0] == "":
            continue

        clust_id = row[0].replace(':','')
        counts = {}
        include = True
        seqs_to_save = []

        for batch in row[1:]:
            seq_names = batch.split(',')
            if len(seq_names) != 1:
                include = False
            if len(seq_names[0]) == 0:
                include = False

            seqs_to_save.extend(seq_names)

        if include:

            print 'saving cluster %s' % clust_id
            out_file = '%s/%s.fasta' % (opts.out_dir, clust_id)
            out_seqs = []

            for seqid in seqs_to_save:
                seq = seqs[seqid.strip()]

                seq.id = os.path.splitext(seq.source)[0]
                out_seqs.append(seq)

            SeqIO.write(out_seqs, out_file, 'fasta')
            file_list.append(out_file)

    print 'Done. Found %s candidate clusters.' % str(len(file_list))
    return file_list

def make_alignments(infiles, prottest, work_dir, numthreads):

    outfiles = []

    for file in infiles:

        print 'working on %s' % file

        out_file = file + '.aligned.fasta'

        out_h = open(out_file, 'w')
        command = 'mafft --quiet --amino --auto --inputorder --thread %s %s' % (numthreads, file)

        output = subprocess.check_output(command.split(' '), shell=False)
        out_h.write(output)
        out_h.close()

        prottest_com = 'java -jar %s -i %s -threads %s -log disabled -JTT -MtREV -DCMut -WAG -RtREV -CpREV -VT -Blosum62 -MtMam -LG -Dayhoff' % (prottest, out_file, numthreads)

        prottest_output = subprocess.check_output(prottest_com.split(' '), shell=False)
        model = ''

        for line in prottest_output.splitlines():

            if line.startswith('Best model according to LnL:'):
                model = line[29:].rstrip()
                found = model.find('+')
                if found != -1:
                    model = model[:found]

        phylip_file = out_file + '.phy'


        convert(out_file, 'fasta', phylip_file, 'phylip-sequential')

        phyml_com = 'phyml --quiet -i %s -m %s --sequential -d aa -b -4 --no_memory_check' %(phylip_file, model)

        subprocess.check_call(phyml_com.split(' '), stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=False)

        outfiles.append(phylip_file+'_phyml_tree.txt')

    return outfiles

def get_alns(in_clade, score_list):

    if in_clade.confidences:
        score = in_clade.confidences[0].value
        score_list.append(score)

    for clade in in_clade.clades:
        get_alns(clade, score_list)

def eval_trees(trees):
    results = []
    for file in trees:

        scores = []
        tree = read(file, 'newick')
        xtree = tree.as_phyloxml()
        get_alns(xtree.clade, scores)
        mean = float(sum(scores)) / float(len(scores))
        minimum = min(scores)

        data = {'mean':mean, 'tree': file, 'min':minimum}
        results.append(data)
    sorted_results = sorted(results, key=lambda m: m['mean'], reverse=True)

    out_h = open('mirlo_report.txt', 'w')
    out_h.write('tree file\tsignal\n')

    for i in sorted_results:
        out_h.write(i['tree']+'\t'+str(i['mean'])+'\n')

if __name__ == '__main__':
    opts = parse_cli()

    try:
        os.stat(opts.out_dir)
    except:
        os.mkdir(opts.out_dir)

    seqs, num_spp = parse_fasta(opts.seq_dir)

    cluster_files = parse_clusters(seqs, opts.clust_file, num_spp)

    num_threads = int(opts.threads)

    print 'starting alignments'

    trees = make_alignments(cluster_files, opts.prottest, opts.out_dir, opts.threads)

    eval_trees(trees)

