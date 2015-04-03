#!/usr/bin/env python
import os
#import pdb
from Bio.Phylo import read
import argparse

def parse_cli():
    parser = argparse.ArgumentParser(description="Calculate phylogenetic signal in a directory of trees")

    parser.add_argument('-i', type=str, dest='work_dir', default='alignments',
                        help='directory containing phylogenetic trees')

    opts = parser.parse_args()
    return opts

def get_alns(in_clade, score_list):

    if in_clade.confidences:
        score = in_clade.confidences[0].value
        score_list.append(score)

    for clade in in_clade.clades:
        get_alns(clade, score_list)


if __name__ == '__main__':
    opts = parse_cli()

    for file in os.listdir(opts.work_dir):
        if file.endswith('_phyml_tree.txt'):
            scores = []
            tree = read(opts.work_dir + file, 'newick')
            xtree = tree.as_phyloxml()
            get_alns(xtree.clade, scores)
            mean = sum(scores) / len(scores)
            minimum = min(scores)
            print "%s\t%s\t%s\t%s" % (file, str(mean), str(minimum), ' '.join(map(str, scores)) )
