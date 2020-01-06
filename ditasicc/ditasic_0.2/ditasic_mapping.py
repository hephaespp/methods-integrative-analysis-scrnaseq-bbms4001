#!/usr/bin/env python

import sys
import os
import glob
import optparse
import time

import numpy as np
from core import ditasic_utils
import json

# Multiprocessing extension
from multiprocessing import Pool


def parallelSimulate(argTup):
    ref,sim,simulator, length, number = argTup
    ditasic_utils.run_simulator[simulator](ref, sim, length, number)
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
    usage = "Usage: %s [OPTIONS] NAMES READS"%(sys.argv[0])
    parser = optparse.OptionParser(usage=usage)
    #parser.add_option('-r', '--reference', type='string', dest='ref', default='./ref/%s.fasta', help='Reference sequence file pattern for the read simulator. Placeholder for the name is "%s". [default: %default]')    
    parser.add_option('-l', '--length', type=int, dest='length', default=100, help='Length of reads used in similarity matrix creation. [default: %default]')
    parser.add_option('-i', '--index', type='string', dest='index', help='Path to kallisto index for all files.') 
    parser.add_option('-t', '--temp', type='string', dest='temp', default='./temp', help='Directory to store temporary simulated datasets and TSV files. [default: %default]')
    # parse arguments
    opt, args = parser.parse_args()
    
    if len(args) != 2: 
        parser.print_help()
        sys.exit(1)

    # create temporary directory if necessary
    if not os.path.exists(opt.temp):
        os.makedirs(opt.temp)

    reads = args[1]
    
    numArgs = len(args)
    threads = 8 

    kal_index = opt.index
    out_path = opt.temp

    contigs = []

    files = [ line.strip() for line in open(args[0]).readlines() ]
    fileNames = [os.path.basename(f) for f in files]

    # Count the contigs in each reference. 
    for i,ref in enumerate(files):
        fai_path = ref+".fai"
        
        # count via fasta index if possible
        if (os.path.exists(fai_path)):
            print("Reading contigs from .fai file %s"%fai_path)
            with open(fai_path) as fai_file:
                #contigFields = [ line.split()[0] for line in fai_file]
                contigs.extend([i]*len(fai_file.readlines()))
        # else just go through the fasta file
        else:
            with open(ref) as f:
                for line in f:
                    if line.startswith(">"):
                        contigs.append(i)

    out_file = out_path + '/' + os.path.splitext(reads)[0].split("/")[-1]

    ditasic_utils.run_kallisto(kal_index, reads, out_file, opt.length)
    vec_unique, vec_shared = ditasic_utils.analyse_tsv(out_file, contigs)

    total = 0
    with open(out_file + "/run_info.json") as json_data:
        data = json.load(json_data)
        total = int(data["n_processed"])               

    print("\nMapping done.\n")
    np.save(os.path.splitext(reads)[0].split("/")[-1]+"_mapped_counts.npy", np.array([x+y for (x,y) in zip(vec_shared, vec_unique)]).astype(np.float64))
    np.save(os.path.splitext(reads)[0].split("/")[-1]+"_total.npy", np.array([total]).astype(np.float64))
