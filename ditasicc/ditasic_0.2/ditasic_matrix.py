#!/usr/bin/env python

import sys
import os
import glob
import optparse
import time

import numpy as np
from core import ditasic_utils

# Multiprocessing extension
from multiprocessing import Pool

def parallelSimulate(argTup):
    ref,sim,ref_probs,length,number = argTup
    ditasic_utils.run_simulator["mason_illumina"](ref, sim, length, ref_probs, number)
    return 1

def parallelMapRow(argTup):
    sam_file, mapper, length, index_file, sim_file = argTup
    param = ""
    if mapper=="kallisto":
        param = "%s"%length
        sam_file = '.'.join(sam_file.split('.')[:-1])

    ditasic_utils.run_kallisto(index_file, sim_file, sam_file, param)
    return 1

if __name__=="__main__":
    usage = """%prog [options] FILES 

    Calculate the similarity matrix.

    First, a set of reads is simulated for every reference genome using the mason
    read simulator.
    Second, the simulated reads of each species are mapped against all reference
    genomes using the Kallisto read mapper.
    Third, the resulting tsv-files are analyzed to calculate the similarity
    matrix. The similarity matrix is stored as a numpy file (-o).

    Input:
    FILES: Files with absolute paths to the reference sequences that we use. 

    See the provided LICENSE file or source code for license information.
    """
    
    parser = optparse.OptionParser(usage=usage)

    parser.add_option('--startprob', type=float, default=0.002, help='Mismatch probability at the first base.')
    parser.add_option('--avgprob', type=float, default=0.004, help='Average per base mismatch probability.')
    parser.add_option('--endprob', type=float, default=0.012, help='Mismatch probability at the last base.')
    
    parser.add_option('--prob-file', type='string', default=None, help='File with mismatch (start, avg, end) probabilities for each reference.')

    parser.add_option('-l', '--length', type=int, dest='length', help='Length of reads used in similarity matrix creation (Required).')
    parser.add_option('-n', type=int, dest='number', default=250000, help='Number of reads sampled from each reference. [default: %default]')
    parser.add_option('-i', '--index', type='string', dest='index', help='Path to existing kallisto index for all files (has to be constructed in same order as the paths in FILES).')    
    parser.add_option('-t', '--temp', type='string', dest='temp', default='./temp', help='Directory to store temporary simulated datasets and TSV files. [default: %default]')
    parser.add_option('--skip', action='store_true', dest='skip', default=False, help='Skip read simulation (assume existing reads in temp) [default: %default].')
    parser.add_option('-o', '--output', type='string', dest='out', default='similarity_matrix.npy', help='Output similarity matrix file. [default: %default]')
    # parse arguments
    opt, args = parser.parse_args()
   
    if len(args) != 1 or opt.length == None: 
        parser.print_help()
        sys.exit(1)

    numArgs = len(args)
    threads = 8 

    kal_index = opt.index
    #reference_path = opt.ref
    out_path = opt.temp

    files = [ line.strip() for line in open(args[0]).readlines() ]
    fileNames = [os.path.basename(f) for f in files]
    #print(fileNames)
    #references = [ reference_path%name for name in files]

    sim_files = [ out_path + "/" + name + ".fastq" for name in fileNames]

    if (opt.prob_file):
        reference_probs = [ [float(y) for y in x.split()] for x in open(opt.prob_file) ]
    else:
        reference_probs = [ [opt.startprob, opt.avgprob, opt.endprob] for name in files ]

    contigs = []

    if kal_index == None:
        print("No index supplied. Creating one ..")
        kal_index = ditasic_utils.build_kallisto_index(files)

    # count number of contigs in the references
    for i,ref in enumerate(files):
        fai_path = ref+".fai"
        if (os.path.exists(fai_path)):
            print("Reading contigs from .fai file %s"%fai_path)
            with open(fai_path) as fai_file:
                #contigFields = [ line.split()[0] for line in fai_file]
                contigs.extend([i]*len(fai_file.readlines()))
        else:
            with open(ref) as f:
                for line in f:
                    if line.startswith(">"):
                        contigs.append(i)

    # create temporary directory if necessary
    if not os.path.exists(opt.temp):
        print("Creating temporary directory %s"%opt.temp)
        os.makedirs(opt.temp)

    n_seq = len(files)
    rng = range(n_seq)

    mapped_reads_unique = np.zeros((n_seq, n_seq))
    mapped_reads_shared = np.zeros((n_seq, n_seq))

    if not opt.skip:
        print("\n--- simulating reads for each reference genome with %s\n"%"Mason_simulator")
        # generate reads for every reference genome
        arglist = [(files[i], sim_files[i], reference_probs[i], opt.length, opt.number) for i in rng]

        if threads==1:
            map(parallelSimulate,arglist)
        else:
            p = Pool(threads)
            p.map_async(parallelSimulate, arglist).get(9999999)
            p.close()

    # map each set of simulated reads onto the index of all sequences
    for i,sim in enumerate(sim_files):
        path = os.path.splitext(sim.split("/")[-1])[0]
        out_file = out_path + '/' + path
        
        # map and output tsv
        ditasic_utils.run_kallisto(kal_index, sim, out_file, opt.length)

        # analyse tsv output
        vec_unique, vec_shared = ditasic_utils.analyse_tsv(out_file, contigs)
        mapped_reads_unique[i, :] = np.array(vec_unique)
        mapped_reads_shared[i, :] = np.array(vec_shared)
    
    np.save(opt.out, mapped_reads_unique + mapped_reads_shared)
