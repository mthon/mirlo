#!/usr/bin/env python
import os
import sys
import pdb
from Bio.AlignIO import convert
from Queue import Queue
from threading import Thread
from time import sleep
import argparse

#work_dir = 'alignments/'

def parse_cli():
    parser = argparse.ArgumentParser(description="Construct phylogentic trees with PhyML")

    parser.add_argument('-i', type=str, dest='work_dir', default='alignments',
                        help='directory containing aligned sequences')

    parser.add_argument('-p', type=str, dest='prottest',
                        help='path to prottest jar file')

    parser.add_argument('-t', type=str, dest='threads',
                        help='number of parallel processes')

    opts = parser.parse_args()

    return opts


def make_tree(q, prottest):

    while not q.empty():
        tree_file = q.get()
        #print 'working on file ' + tree_file
        prottest_com = 'java -jar %s -i %s/%s  -threads 1' % (opts.prottest, opts.work_dir, tree_file)
        # -JTT -MtREV -DCMut -WAG -RtREV -CpREV -VT -Blosum62 -MtMam -LG -DayHoff
        # -all-matrices -all-distributions
        #print 'running: ' + prottest_com
        output = os.popen(prottest_com)

        model = ''
        for line in output:
            if line.startswith('Best model according to LnL:'):
                model = line[29:].rstrip()
                found = model.find('+')
                if found != -1:
                    model = model[:found]
        phylip_file = opts.work_dir + '/' + tree_file+'.phy'
        #pdb.set_trace()
        print phylip_file + '\t' + model
        convert(opts.work_dir + '/' + tree_file, 'fasta', phylip_file, 'phylip-sequential')

        phyml_com = 'phyml -i %s -m %s --sequential -d aa -b -4 --no_memory_check' %(phylip_file, model)
        print 'running: ' + phyml_com
        #pdb.set_trace()

        os.popen(phyml_com)
        print 'phyml is done'
    q.task_done()


if __name__ == '__main__':
    opts = parse_cli()

    tree_queue = Queue()
    workers = []

    for file in os.listdir(opts.work_dir):
        if file.endswith('aligned.fasta'):
            tree_queue.put(file)

    num_threads = int(opts.threads)

    for i in range(num_threads):
        worker = Thread(target=make_tree, args=(tree_queue,opts))
        worker.setDaemon(True)
        worker.start()

    while not tree_queue.empty():
        sleep(1)
