#!/usr/bin/env python
import os
import sys
import shutil
import pdb
import subprocess
import argparse
from Bio import SeqIO
from Queue import Queue
from threading import Thread
import threading
from time import sleep

def parse_cli():
    parser = argparse.ArgumentParser(description="split a fasta file into segments and run a self-blast")
    parser.add_argument('-i', required=True, type=str, dest='infile', help='input fasta file')
    parser.add_argument('-t', type=int, dest='threads', help='number of parallel threads')
    parser.add_argument('-e', type=str, dest='evalue', help='e-value threshold')
    parser.add_argument('-o', type=str, dest='outfile', help='output file')
    parser.add_argument('-w', type=str, dest='work_dir', help='temporary working directory', default='blast/')

    args = parser.parse_args()
    return args

def split_fasta(file, num_threads, work_dir):
    seq_dict = SeqIO.index(file, 'fasta')
    num_seqs = len(seq_dict)
    seqs_per_thread = num_seqs / num_threads+1
    #pdb.set_trace()
    counter = 0
    seq_batch_num = 0
    seqs_batch = []

    file_h = open(file, 'rU')
    batch_file_names = []
    for seq in SeqIO.parse(file_h, 'fasta'):
        #if counter <= seqs_per_thread:
        seqs_batch.append(seq)
        counter += 1
        if counter== seqs_per_thread:
            out_file_name = work_dir+'/batch-'+str(seq_batch_num)+'.fasta'
            SeqIO.write(seqs_batch, out_file_name, 'fasta')
            seqs_batch = []
            seq_batch_num +=1
            counter = 0
            batch_file_names.append(out_file_name)

    # write the last set of seqs to a file
    out_file_name = work_dir+'/batch-'+str(seq_batch_num)+'.fasta'
    SeqIO.write(seqs_batch, out_file_name , 'fasta')
    batch_file_names.append(out_file_name)
    return batch_file_names


def make_blast_db(file):
    # format the blast database
    command = ['makeblastdb', '-in', file, '-dbtype', 'prot']
    subprocess.check_call(command)

def run_blast(split_files, blast_db, work_dir, evalue, outfile):

    def do_it(queue, command):
        one_year = 365 * 24 * 60 * 60
        file = queue.get(True, one_year)
        print 'running blast on file '+file+'\n'
        print ' '.join(command) + '\n'
        subprocess.check_call(command)
        print 'blast done: ' + file
        queue.task_done()

    blast_queue = Queue()
    workers = []
    out_files = []

    for file in split_files:
        blast_queue.put(file)

        split_outfile = file+'.blast'
        blast_command = ['blastp','-query',file,'-db', blast_db, '-evalue', evalue, '-out', split_outfile, '-outfmt', '7']
        worker = Thread(target=do_it, args=(blast_queue,blast_command,))
        worker.setDaemon(True)
        worker.start()
        # pdb.set_trace()
        out_files.append(split_outfile)

#     while not blast_queue.empty():
#         sleep(1)
#     print 'queue size: ' + str(blast_queue.qsize())
    #blast_queue.join()

    isRunning = True
    while isRunning:
        sleep(1)
        thread_count = len(threading.enumerate())

        print "running on "+str(thread_count)+" threads"
        if thread_count ==1:
            isRunning = False

    outfile_h = open(outfile, 'w')
    for this_file in out_files:
        this_file_h = open(this_file)
        for line in this_file_h:
            outfile_h.write(line)
            #pdb.set_trace()

if __name__ == '__main__':
    cli = parse_cli()

    infile_basename = os.path.basename(cli.infile)
    infile_workdir = os.path.join(cli.work_dir, infile_basename)

    # check if work dir exists
    if not os.path.isdir(cli.work_dir):
        os.mkdir(cli.work_dir)

    shutil.copy(cli.infile, infile_workdir)

    split_files = split_fasta(cli.infile, cli.threads, cli.work_dir)

    make_blast_db(infile_workdir)

    run_blast(split_files, infile_workdir, cli.work_dir, cli.evalue, cli.outfile)

