#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## ------------------------------------------------------------------------------------------------------------------------------------ ## 
###############################################################################
##     COUNTHIC - Counts reads within txt files made by HicLite in SLURPY    
###############################################################################
## ------------------------------------------------------------------------------------------------------------------------------------ ## 
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ------------------------------------------------------------------------------------------------------------------------------------ ## 
##      MODULE LOADING
## Bring in pandas 
import pandas as pd, sys, json, gzip
## Load in defaults 
from defaults import sortglob, aligndir, dictzip, diagdir, isgzip
## Loadin the correct backend for matplotlib
import matplotlib
## Set the needed backend 
matplotlib.use('Agg')
## Load in matplot lib 
from matplotlib import pyplot as plt
## Bringin ftns from os
from os.path import getsize, exists
## ------------------------------------------------------------------------------------------------------------------------------------ ## 

## ------------------------------------------------------------------------------------------------------------------------------------ ## 
##      FUNCTION SETTING
def countfile(inpath:str) -> int: 
    """Counts the lines within file."""
    ## Open the file for reading
    with open(inpath,'rb') as infile:
        ## Iterate thru the handle 
        for count,line in enumerate(infile):
            pass 
    ## Return the count + 1
    return count + 1

def countgzip(inpath:str) -> int:
    """Counts the lines within a gzipped file."""
    with gzip.open(inpath, 'rb') as infile:
        ## Iterate over the file
        for count, line in enumerate(infile):
            pass 
    return count + 1

## Ftn for summing rows across chunks 
def linecount(inpath:str) -> int:
    """Returns the number of rows within a file (can be gzipped) or returns zero if file is empty."""
    ## Return the sum of the number of rows in the file or zero if empty 
    return (countgzip(inpath) if isgzip(inpath) else countfile(inpath)) if (getsize(inpath) and exists(inpath)) else 0
    
## Ftn for parsing json file
def parsejson(injson:str) -> list:
    """Load and parses an input json file, returns object as json dictionary."""
    ## Loads the json file
    with open(injson,'r') as infile:
        data = json.load(infile)
    ## Return the parsed data
    return data

## Ftn for getting total read counts
def jsontotals(inpath:str) -> int:
    """From a .json file produced by fastp, gathers the total read counts before filtering."""
    ## Return the total
    return int(parsejson(inpath)['summary']['before_filtering']['total_reads'])/2
## ------------------------------------------------------------------------------------------------------------------------------------ ## 

## ------------------------------------------------------------------------------------------------------------------------------------ ## 
##      VARIABLE SETTING
## Set colors used in plotting 
colors_to_use =  ['tab:blue','tab:grey','tab:green','tab:orange','black','tab:purple','tab:brown','tab:red']
## Bring in cycle 
from itertools import cycle       
## Make a color
mycolors = cycle(colors_to_use)
## Set the fontsize 
myfs = 12
## Set the line width
linewidth = 30
## ------------------------------------------------------------------------------------------------------------------------------------ ## 

## ------------------------------------------------------------------------------------------------------------------------------------ ## 
##      BODY of MAIN EXECUTABLE 
## If the script is called 
if __name__ == "__main__":
    ## Set the run name 
    run_name = sys.argv[1]

    ## Gather the txt files 
    txtfiles = sortglob(f'./{aligndir}/*.txt') + sortglob(f'./{aligndir}/*.txt.gz')
    ## Check our work
    assert len(txtfiles), "ERROR: No input txt or txt.gz files were detected"

    ## Gather the json files 
    jsonfiles = sortglob(f'./{diagdir}/*.json')

    ## Calc the totals 
    totals = [jsontotals(p) for p in jsonfiles]

    ## Gather the jason short name
    json_short_names = [j.split('/')[-1].split('.fastp.')[0] for j in jsonfiles]

    ## Gather the basename of the files 
    basefiles = [t.split('/')[-1] for t in txtfiles]

    ## Initilze list counts per txt file, add the totals per json file, and the summed totals 
    counts_per_file = [linecount(inpath) for inpath in txtfiles] + totals + [sum(totals)]

    ## Get the short names
    short_names = ['.'.join(f.split('.txt')[0].split('.')[1:]) for f in basefiles] + json_short_names + ['Total']

    ## Make a dicitonary of counts and colors by name
    name_counts = dictzip(short_names,counts_per_file)
    name_colors = dictzip(short_names,mycolors)

    ## Call a figure set face color to wight 
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    fig.set_facecolor('w')

    ## Plot the lines of each count 
    [plt.vlines(i, 0, name_counts[sname], label=sname,color = name_colors[sname], linewidth=linewidth) 
                for i,sname in enumerate(short_names)]
    ## Add x ticks
    plt.xticks(range(len(short_names)), short_names, fontsize=myfs, rotation=90)
    ## Add x and y label
    plt.xlabel('Hi-C File', fontsize=myfs)
    plt.ylabel('Number of Fragments\n', fontsize=myfs)
    ## Save figure 
    plt.savefig(f'./{diagdir}/{run_name}.hic.counts.png',dpi=300,bbox_inches='tight')

    ## Make the counts into a csv 
    dfout = pd.DataFrame(zip(basefiles + json_short_names + ['Total'],short_names,counts_per_file),
                         columns=['File','Shortname','Fragments'])
    ## Saveout the df
    dfout.to_csv(f'./{diagdir}/{run_name}.hic.counts.csv',index=False)
## ------------------------------------------------------------------------------------------------------------------------------------ ## 
## End of file 