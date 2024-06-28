#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
###################################################
##  READ SPLITER - splits a bam file (v 3.0.0)   ##
###################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ---------------------------------------------- LOAD IN MODULES ------------------------------------------------- ## 
## Load in ftns from pysam tools 
from pysamtools import loadbam, checksam, writeset, ifprint, outnames

## Bring in pandas dataframe
import pandas as pd 

## ------------------------------------------- SET DEFUALT VARIABLES ---------------------------------------------- ##
## Set version
version = '3.0.0'

## Set the description
split_desc = "READ SPLITER (v %s): Splits an input bam file from a WGS, ChIP, ATAC experiment on mapped, placed, unmapped, and mtDNA read sets."%version

## Define the error messages
write_read_err  = "ERROR: Unable to write %s read names to file!"

## Set warning messages 
read_sort_warn  = "WARNING: There was an unknown error in sorting the reads; likely due to non-unique read names."

## Set the mapping quality and mito contig
mito_contig = 'chrM'

## Set help messgaes
B_help = "Path to input bam file from Hi-C experiment to split on genomic regions."
M_help = "Name of the mitochondrial contig (default: %s)."%mito_contig

## ------------------------------------------- DEFINE FUNCTIONS ------------------------------------------------ ##
## Ftn for making table 
def recordstable(pyobj, colnames = ['Readname','Isread1','Unmapped','Chromosome','Mapq']):
    """Generates a table of basic information from pyrces, including query name, is read one, and if pairs are unmapped."""
    ## Returns a dataframe with mapping info per record
    return pd.DataFrame([(r.query_name,r.is_read1,r.is_unmapped,r.reference_name,int(r.mapq)) for r in pyobj] , columns=colnames)

## Ftn returns set of all reads
def readset(df):
    """Return the set of Readnames"""
    ## Return the set of readnaems
    return set(df.Readname.tolist())

## ----------------------------------------- BODY of EXECUTABLE ----------------------------------------------- ## 
## if the script is envoked
if __name__ == "__main__":
    ## Check the versions of samtools
    assert checksam(), 'ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again.\n'

    ## Load argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = split_desc)

    ## Add the required arguments and default arguments 
    parser.add_argument("-b", "--bam", "--input", dest="B", type=str, required=True,  help=B_help, metavar='./path/to/bam')
    parser.add_argument("-M", "--mitochondria",   dest="M", type=str, required=False, help=M_help, metavar=mito_contig, default=mito_contig)

    ## set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set input variables, input bam file path, mito contig name, and mapping quality threshold 
    inbampath, mito = inputs.B, inputs.M

    ## Format the name of output files generated by this script 
    mapped_name, placed_name, mito_name, unmapped_name, bedpe_file = outnames(inbampath,mito)

    ## Bring in remove ftn 
    from os import remove

    ## Load in is file from os
    from os.path import isfile

    ## Remove old runs if needed
    [(remove(pr) if isfile(pr) else None) for pr in [mapped_name, placed_name, mito_name, unmapped_name]]
       
    ## Load in py recs and make a data tabel 
    rectable = recordstable(loadbam(inbampath))

    ## Make a set of all read names
    allreads = readset(rectable)

    ## Parse the mitochondrial reads if mito was set 
    mitocondrial = readset(rectable[(rectable.Chromosome==mito)]) if mito else set()

    ## Write the mitochondrial name to file 
    if mito: ## Assert the file is made 
        assert isfile(writeset(mito_name, mitocondrial)), write_read_err%mito

    ## Gather read1 and read2 dataframes without mito mapping reads
    read1, read2 = rectable[(rectable.Isread1) & ~(rectable.Readname.isin(mitocondrial))], rectable[~(rectable.Isread1) & ~(rectable.Readname.isin(mitocondrial))]

    ## Gather sets of unmapped r1 and r2
    unmapped_r1, unmapped_r2  = readset(read1[(read1.Unmapped)]), readset(read2[(read2.Unmapped)])

    ## Gather intersection were both reads unmapped
    unmapped = unmapped_r1 & unmapped_r2

    ## Gather the unmapped reads another way to check our way
    unmapped_2 = readset(rectable[(rectable.Chromosome.isnull())])

    ## Check our work
    ifprint('WARNING: The number of unmapped reads did not match during parsing.', len(unmapped - unmapped_2) == 0)

    ## Write the unmapped read names to file 
    assert isfile(writeset(unmapped_name, unmapped)), write_read_err%'unmapped'

    ## Gather placed reads
    placed = (unmapped_r1 | unmapped_r2) - unmapped

    ## Write the placed read names to file 
    assert isfile(writeset(placed_name, placed)), write_read_err%'placed'

    ## Segregate the reads to parse 
    mapped = allreads - (unmapped | placed | mitocondrial)

    ## Write the mapped read names to file
    assert isfile(writeset(mapped_name, mapped)), write_read_err%'mapped'

    ## Group the sets
    all_sets = [mitocondrial, unmapped, placed, mapped]

    ## Check our owrk, iterate thru the sets
    for i, s1 in enumerate(all_sets):
        for j, s2 in enumerate(all_sets):
            if j>i: ## Make sure the intersection of the sets is zero!
                ifprint(read_sort_warn, (len(s1.intersection(s2)) != 0))
## End of file 