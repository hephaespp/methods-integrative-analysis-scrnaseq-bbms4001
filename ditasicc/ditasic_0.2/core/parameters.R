# Copyright (c) 2016,
# Martina Fischer, fischerm@rki.de, Robert Koch Institute, Germany, 
# All rights reserved. For details, please note the license.txt.

# Function for parameter parsing and preparation, 
# package loading and error handling


loadPackages <- function(names, quietly=TRUE){
  # Tries to require R packages/libraries and stop executions if failing
  #
  # Args:
  #   names: vector with package names
  # Return:
  #   Nothing, but will stop execution if packages are missing
  
  check <- vector(mode="logical", length=length(names))
  
  for (i in 1:length(names)){
    suppressMessages(check[i] <- require(names[i], 
                                         quietly=quietly, 
                                         character.only=TRUE))
  }
  
  if (!all(check)){
    stop("missing R-package(s)", names[!check])
  }
}

# DiTASiC argument parser

#command args:  ref_path, mat_path, counts_s1_path, N1_path, counts_s2_path, N2_path, 
#optional:      filtering=T/F, pval.threshold, min.threshold 

parseParameters <- function(args=commandArgs(trailingOnly = TRUE)){
  # Parses parameters from command line to DiTASiC.RScript
  # 
  # Args:
  #   args: commandline parameters, default = commandArgs(trailingOnly = TRUE) 
  # Returns:
  #   List of parameter values as provided by parse_args (optparse package)
  
  parser.options <- list(
    
    # Input parameters required
    make_option(opt_str=c("-r", "--refs"),
                action="store",
                help="Indicate the taxa names file containing the 
                abolute path to all considered taxa references in this analysis 
                (already required in previous 'ditasic_mapping' and 'ditasic_matrix')", type = 'character', default=NULL),
    make_option(opt_str=c("-a", "--Mat"),
                action="store",
                help="Indicate path to the similarity_matrix.npy file (output of 'ditasic_matrix').", type = 'character', default=NULL),
    make_option(opt_str=c("-x", "--counts_s1"),
                action="store",
                help="Indicate path to the mapped count vector of sample 1: 
                sample.npy file (default output of 'ditasic_mapping').", type = 'character', default=NULL),
    make_option(opt_str=c("-y", "--counts_s2"),
                action="store",
                help="Indicate path to the mapped count vector of sample 2: 
                sample.npy file (default output of 'ditasic_mapping').", type = 'character', default=NA),
    make_option(opt_str=c("-n", "--N_s1"),
                action="store",
                help="Indicate path to the total vector of sample 1, containing the total number
                of input reads: total.npy file (default output of 'ditasic_mapping').", type = 'character', default=NULL),
    make_option(opt_str=c("-m", "--N_s2"),
                action="store",
                help="Indicate path to the total vector of sample 2, containing the total number
                of input reads: total.npy file (default output of 'ditasic_mapping').", type = 'character', default=NA),
    
    # optional parameters
    make_option(opt_str=c("-f", "--filter"), 
                action="store",
                default=F,
                type="logical",
                help="Indicate whether a filtering should be applied to detect and 
                remove false-positive taxa in the data.(default = F)"),
    make_option(opt_str=c("-o", "--output"), 
                action="store",
                default="DiffAbund_Result.txt",
                type="character",
                help="Indicate the name of the output file"),
    make_option(opt_str=c("-p", "--pval_thres"), 
                action="store",
                default=0.05,
                type="numeric",
                help="Indicate a p-value threshold to remove false-positive taxa
                in case of filtering applied.(default = 0.05)"),
    make_option(opt_str=c("-t", "--min_thres"), 
                action="store",
                default=0,
                type="integer",
                help="Indicate a minimum number of reads to assign significant 
                taxa existence.(default = 0)"),
    make_option(opt_str=c("-s", "--seed"), 
                action="store",
                default=1448,
                type="integer",
                help="Indicate a seed for sampling processes used in creating the 
                empirical distributions.(default = 1448)")
  )
  
  my.parser <- OptionParser(usage=
                             "usage: ditasic -r REFPATH_FILE -a SIMILARITY_MATRIX -x COUNTS_S1 -n TOTAL_S1 [-y COUNTS_S2 -m TOTAL_S2] [options]",
                           option_list=parser.options,
                           prog = NULL, description = "", epilogue = "")
  
  parameters <- parse_args(object = my.parser,
                           args = args,
                           print_help_and_exit = TRUE,
                           positional_arguments = FALSE)
  
  return(parameters)
  
}
# end
