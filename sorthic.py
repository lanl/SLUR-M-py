#!/usr/bin/env python
#################################################
##      Hi-C Chromosome Sorter  (v 0.0.0)      ##
#################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Set version
version = '0.0.0'

## Set the description
sort_desc = "Hi-C Chromosome sorter (v %s): Sorts an input txt file (space delimenated) representing Hi-C contacts."%version

## Set help messages 
i_help = "Paths to Hi-C text files for sorting."
o_help = "Path to output save file."

## Load in short cols and dict 
from splithic import short_columns, short_type_dict

## Bring in pandas
import shutil

## Bring in exists form os 
from os.path import exists, getsize

## ------------------------------------------------------ BODY of MAIN EXECUTABLE --------------------------------------------------------- ##
## If the script is called 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = sort_desc)
    ## Add the required arguments
    parser.add_argument("-i", dest="I", nargs='+', required=True, help=i_help)
    parser.add_argument("-o", dest="O", type=str,  required=True, help=o_help)
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Get input paths assuming they are sorted
    inputpaths = inputs.I
    outputpath = inputs.O 

    ## Set the ouput duplicates file 
    outdups = outputpath.split('.tohic')[0] + '.duplicates.txt'

    ## Set the columns 
    used_cols = short_columns[:8]
    used_dict = [short_type_dict[c] for c in used_cols]

    ## Initlize listof duplicate paths 
    dups_paths = []

    ## Open the output file for writing
    with open(outputpath,'bw') as outhicfile:
        ## Itreate over the path
        for i,inpath in enumerate(inputpaths):
            
            ## Gather the duplicates for this file, set the dup path 
            chrom_dups = inpath.split('.txt')[0] + '.duplicates.txt'
            ## Append dups pahs
            dups_paths.append(chrom_dups)

            ## If the file exists and has a file size greater than zero 
            if exists(inpath) and getsize(inpath):
                ## Open for reading 
                with open(inpath,'rb') as inhicfile:
                    ## Copy out to out hic file 
                    shutil.copyfileobj(inhicfile,outhicfile)

    ## Open the output duplicates file for writing 
    with open(outdups,'bw') as destination:
        ## Iterate over the path
        for i,inpath in enumerate(dups_paths):
            ## If the file exists and has a file size greater than zero 
            if exists(inpath) and getsize(inpath):
                ## Open duplicate for reading 
                with open(inpath,'rb') as source:
                    ## Copy to file 
                    shutil.copyfileobj(source, destination)
## End of file