#! /usr/bin/env python3

import os
import sys
import random
import subprocess

def analyse_tsv(tsv_dir, contig_index):
    ref_size = contig_index[-1]+1
    uniq_n = len(contig_index)
    
    vec_unique = [0.]*ref_size
    vec_shared = [0.]*ref_size

    # pseudoalignments mapped to sets of contigs
    with open(tsv_dir + "/pseudoalignments.tsv") as tsv:
        pseudo = [ int(line.split()[1]) for line in tsv.readlines() ]

    # equivalence classes 
    ec_index = []
    with open(tsv_dir + "/pseudoalignments.ec") as ec:
        ec_index = [ (line.split()[1]).split(",") for line in ec.readlines() ]

    # sum the equivalence classes as they distribute reads among the reference contigs
    for i, ec in enumerate(ec_index):
        # the contig received this many hits
        contig_count = pseudo[i]
        
        # for each of the contigs covered by the EC got the hits
        processed_refs = set()
        added = []
        for contig in ec:
            # the hits go to this reference
            ref = contig_index[int(contig)]

            # only add to each reference once
            if ref in processed_refs:
                continue

            processed_refs.add(ref)
            added.append( (ref, contig_count) )

            if i < uniq_n:
                vec_unique[ref] += contig_count
            else:
                vec_shared[ref] += contig_count

        assert( len(processed_refs) == len(added) )

    return(vec_unique, vec_shared)

def build_kallisto_index(files):
    # NOTE: let user name index when newly created?
    index_name = "kalindex_%iref"%(len(files))
    while os.path.exists(index_name): # do not overwrite existing indices
        raise Exception("Error in automatic index creation: Index '%s' exists already."%index_name)

    print("Saving index as %s"%index_name)
    file_args = ' '.join("'{0}'".format(f) for f in files)
    command   = ("kallisto index -i %s %s"%(index_name,file_args))

    subprocess.check_call(command, stderr=subprocess.STDOUT, shell=True)

    return index_name

def run_kallisto(index, reads, out, param="100"):
    #if (out.endswith(".bam")):
    #    command = "kallisto pseudo -t 1 -i {index} -o /tmp/kallisto --single -l {param} -s 20 {reads} --pseudobam | samtools view -bS - > {samfile}".format(reads=reads, index=index, samfile=out, param=param)
    #elif (out.endswith(".sam")):
    #    command = "kallisto pseudo -t 1 -i {index} -o /tmp/kallisto --single -l {param} -s 20 {reads} --pseudobam > {samfile}".format(reads=reads, index=index, samfile=out, param=param)
    #else:
    command = "kallisto pseudo -t 32 -i {index} -o '{samfile}' --single -l {param} -s 20 '{reads}'".format(reads=reads, index=index, samfile=out, param=param)

    with open(os.devnull, 'wb') as devnull:
            subprocess.check_call(command, stderr=subprocess.STDOUT, shell=True)
    return 1

# Read Simulators
#
# Define caller functions for the read simulators here. These fuctions call
# the read simulators on the command line with the correspondig parameters.
# Output is a SAM file at the specified location. Each function takes 3
# mandatory string arguments:
#   ref: name of the reference sequence file
#   out: name of the output reads file
#
# Note: since simulators tend to have many tuning parameters, we encourage
# to create a seperate caller function for every scenario.

def run_mason_illumina(ref, out, length, ref_probs, num):
    command = "mason_simulator -n %i --illumina-read-length %i --illumina-prob-mismatch %f --illumina-prob-mismatch-begin %f --illumina-prob-mismatch-end %f -ir '%s' -o '%s'"%(num, length, ref_probs[1], ref_probs[0], ref_probs[2], ref, out)
    
    print("Executing:",command)
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(command, stderr=subprocess.STDOUT, shell=True)
        #subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT, shell=True)

    # remove the needless SAM file
    #command = "rm %s.sam"%out
    #os.system(command)
    return 1

def run_dwgsim(ref, out, length, num):
    command = "dwgsim -c 2 -1 %i -2 0 -r 0 -y 0 -e 0.002 -N %i -f TACG %s %s"%(length, number, ref, out)
    print("Executing:",command)
    os.system(command)
    # remove all additional files and rename reads file
    command = "mv {out}.bfast.fastq {out} && rm {out}.bwa.read1.fastq && rm {out}.bwa.read2.fastq && rm {out}.mutations.txt".format(out=out)
    os.system(command)
    return 1


run_simulator = dict( mason_illumina=run_mason_illumina,
                      dwgsim=run_dwgsim,)
                    
"""
How to add your custom read simulator

1. Create a caller function
   - Create a copy of one of the existing functions, e.g. run_dwgsim
   - Rename it, customize it, but DO NOT TOUCH the interface!

2. Add the caller function to the run_simulator dict
   - The dict entry should have the format: [name] = [caller function]

3. Now you can use your simulator:
 >>> import tools_lib
 >>> tools_lib.run_simulator["name"]
"""
