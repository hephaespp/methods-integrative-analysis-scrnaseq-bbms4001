#!/bin/bash

### An example run of a typical Ditasic use case

# Generate the similarity matrix of the 35 reference files.
# Without the -i parameter given, it will automatically create
# an index for Kallisto with the name "kalindex35_ref". 

ditasic_matrix.py -l 100 -o output/similarity_matrix35.npy data/reference_paths

# if the index exists already - pass with parameter -i 
#./ditasic/ditasic_matrix.py -l 100 -i kalindex_35ref -o output/similarity_matrix35.npy data/reference_paths

# Map the samples to the reference files for the raw abundance values. 
ditasic_mapping.py -l 100 -i kalindex_35ref data/reference_paths data/sample1.fastq
ditasic_mapping.py -l 100 -i kalindex_35ref data/reference_paths data/sample2.fastq

mv sample1_mapped_counts.npy sample1_total.npy output
mv sample2_mapped_counts.npy sample2_total.npy output

# Estimate abundance for sample 1
ditasic -r data/reference_paths -a output/similarity_matrix35.npy \
                  -x output/sample1_mapped_counts.npy -n output/sample1_total.npy \
                  -f F -t 750 -o output/sample1_abundance.txt

# Estimate differential abundance between sample 1 and sample 2
# (New parameters: Mapping counts for two samples as the input)
ditasic -r data/reference_paths -a output/similarity_matrix35.npy \
                  -x output/sample1_mapped_counts.npy -n output/sample1_total.npy \
                  -y output/sample2_mapped_counts.npy -m output/sample2_total.npy \
                  -f F -t 750 -o output/sample1_sample2_diff_abundance.txt
