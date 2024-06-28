#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
filterdesc = "Filters an input SAM file for genomically mapped, unmapped (or placed), and mitochondiral mapping reads."
## ----------------------------- v 0.0.0 ------------------------------ ## 
##   Pre (Hi-C) Filter   
## -------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ----------------------------------- MODULE LOADING ------------------------------------ ##
## Bring in column names of sam file 
from pysamtools import issam, writeset, loadsam

## Bring in ftns and options from slurpy
from defaults import mito, M_help, pathexists, chunks, Z_help, reset

## ------------------------ FUNCTION and VARIABLE FILTERING ------------------------------ ## 
## Set help strings
S_help = "Path to input .sam file from bwa."

## ---------------------------- MAIN SCRIPT & ARGUMENT PARSING -------------------------- ## 
## If the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = filterdesc)
    ## Add the required arguments
    parser.add_argument("-s", dest="S", type=str,  required=True,  help=S_help, metavar='./path/to/input.sam'  )
    ## Add optional variables 
    parser.add_argument("-M", dest="M", type=str,  required=False, help=M_help, metavar='chrM', default=mito   )
    parser.add_argument("-Z", dest="Z", type=int,  required=False, help=Z_help, metavar='n',    default=chunks )
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Patch input vars 
    insampath = inputs.S   ## Set input path
    mito      = inputs.M   ## Set the mito dna contig
    chunksize = inputs.Z   ## Set chunk size for parsing sam file 

    ## Assert we have a sam file
    assert issam(insampath)

    ## Set the output paths 
    path_head = insampath.split('.sam')[0]    ## Gather the head path 
    mito_path = path_head + f'.{mito}.txt'    ## Set the mtDNA mapping reads 
    unmp_path = path_head + f'.unmapped.txt'  ## Set path to the unmapped reads
    save_path = path_head + f'.genomic.txt'   ## Set name of output path to save

    ## Remove previous runs
    reset([mito_path,unmp_path,save_path])

    ## initiate sets of un mapped and mitochondiral mapped reads 
    unmp_reads = set()
    mito_reads = set()

    ## Open an instance with pandas, parse thru once to get the unmapped and mtDNA mapping reads 
    with loadsam(insampath,chunksize) as samchunks:
        ## Iterate the sam chunks 
        for thesam in samchunks:
            ## Parse the chunks to get the set of un mapped reads 
            unmp_reads = unmp_reads | set(thesam[(thesam.Rname=='*') | (thesam.Rnext=='*') | (thesam.Cigar=='*')].Qname)
            ## Gather the set of mito mapping reads 
            mito_reads = mito_reads | set(thesam[(thesam.Rname==mito) | (thesam.Rnext==mito)].Qname)
    
    ## Write the mito mapping reads to file
    writeset(mito_path,mito_reads) if len(mito_reads) else None 
    ## Write the un-mapped reads set to file
    writeset(unmp_path,unmp_reads) if len(unmp_reads) else None 
    
    ## Open an instance with pandas, parse thru a second time to save out the mapped reads 
    with loadsam(insampath,chunksize) as samchunks:
        ## Iterate the sam chunks 
        for thesam in samchunks:
            ## Filter those names from above and remove from data
            savesam = thesam[~((thesam.Qname.isin(unmp_reads)) | (thesam.Qname.isin(mito_reads)))]
            ## Save the filtered, smaller dataframe
            savesam.to_csv(save_path, sep='\t', index=False, header=(not pathexists(save_path)), mode='a' if pathexists(save_path) else 'w') if savesam.shape[0] else None 
## End of file 