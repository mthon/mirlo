Mirlo
=====

Overview
--------

Mirlo is protocol and a set of scripts to help automate the construction
of data sets for  multi-gene phylogenetic analyses. It takes whole
proteomes as input, identifies single-copy gene families, aligns the
proteins, and evaluates the 'phylogenetic signal' of each alignment.

What mirlo does:

* Cluster proteins using BLAST and MCL.

* Identify single-copy gene families (clusters.)

* Align each single-copy family using MAFFT

* Construct a phylogenetic tree of each family using PhyML and compute
SH-like support values for each branch.

* Compute a 'phylogenetic signal' for each family by computing the mean
SH-like values for all branches in each tree. This step is inspired by
Salichos L, and Rokas A: Inferring ancient divergences requires genes
with strong phylogenetic signals. Nature 2013, 497:327–331. The authors
show that genes with higher phylogenetic signals have phylogenies that
are more congruent with the species tree. Here, I use SH-like support
values instead of bootstrap values because they are much faster to
compute.

What mirlo does NOT do:

* It does not edit the alignments

* It does not construct a phylogenetic tree from the concatenated
alignment.

Mirlo is a work in progress.

Contact
-------

If you have questions or need help, email me:

Michael Thon mthon@usal.es


Software Requirements
---------------------

- MCL 
- ProtTest 3.4 https://code.google.com/p/prottest3/
- PhyML
- mafft
- BioPython

All of the above software packages are available on most linux
distributions except for ProtTest.

Installation
------------

Uncompress the Mirlo distribution file somewhere on your computer.

Instructions
------------

You should have one fasta format file of proteins for each species that
you plan to include in your analysis. Put all of the fasta files into
the the same directory.

1. Prepare the file species_list.csv. Column 1 should contain a three
letter code for each species. Column 2 contains the genus and species of
each taxon in your analysis.  Column 3 contains the file name of the
fasta file that you put in the 'fasta' directory.

1. run rename.py.  This script reads the proteins in each file and
renames them, appending the three letter code specified in
species_list.csv The output of this script is a file containing the
renamed proteins from all of the input proteims.

1. run an 'all vs. all' blast of the proteins.

    `makeblastdb -in all.fasta -dbtype prot`

    `blastp  -query all.fasta -db all.fasta -outfmt 7 -num_threads 24
    -evalue 1e-3 -out selfblast.txt`

    * A compute cluster really helps out with this task. If all you have
    is a multi-core server, check out the run_blast.py script which is a
    little more efficient at parallelization than the ncbi blast
    programs. Alternatively, you can use blat, which is much faster, but
    probably less sensitive:

        `blat all.fasta all.fasta -out=blast9 -prot out.blat`

1. Cluster the sequences with MCL

    `mcxdeblast --m9 --line-mode=abc -out=- selfblast.txt | mcl - --abc
    -I 1.2 -o all.clusters`

1. run make_alignments.py This script identifies which clusters contain
one protein from each species and aligns the sequences using mafft.

1. run make_trees.py This script runs prottest to select the best model,
then runs PhyML using the best model. PhyML also computes SH-like
support values for each branch.

1. run evaluate_trees.py This step computes the average SH-like support
value for each tree and prints to stdout. the output is tab separated
and contains the following: the tree file name, the average SH-like
support value, the minimum support value, a space separated list of all
the support values in the tree. This last column is mainly for debugging
purposes.

1. run rename_alignments.py
