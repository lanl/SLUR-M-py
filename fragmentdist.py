#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
###########################################################
##    Fragment Distribution Ploting Library & Script     ##
###########################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ---------------------------------------------- LOAD IN MODULES --------------------------------------------- ## 
## Bring in mods
import numpy as np, pandas as pd 

## Bring in matplot lib
from matplotlib import pyplot as plt 

## Load in ftns from slurpy
from defaults import basenobam

## Load in ftn from pysamtools
from pysamtools import loadbam, isbam, hasix

## -------------------------------------------- GENERAL FUNCTIONS --------------------------------------------- ## 
## Ftn for calculating index sizes 
def indexsizes(inbam):
    """Gathers the index size from input bam files."""
    ## Return the template lengths from int input bam file 
    return np.array([r.template_length for r in inbam if (r.is_read1 and (not r.is_supplementary) and (not r.is_secondary))])

def indexfreq(indxsizes,delta=10,maxsize=1000):
    """Calculates the frequency of the insert sizes."""
    ## Gather the unique counts, this version is old 
    k = np.abs(indxsizes)  #l,p = np.unique(np.abs(insizes),return_counts=True)
    ## Generate the bins
    bins = np.arange(0,maxsize+delta,delta)
    ## Return the unique, binned fragmentsizes and the counts dividied by the total
    return bins, np.array([np.sum((k<b)&(k>=(b-delta))) for b in bins])/len(indxsizes)

## Ftn for plotting fragments
def plotfragments(inpaths,outsavepath=None,maxindexsize=1000,binsize=10,ax=None):
    """Plots fragment length distributions from input bam files"""
    ## Set the axis if it was give
    plt.sca(ax) if ax else None 
    ## Initilize the list of dfs
    dfs = []
    ## Iterate thru the paths 
    for inpath in inpaths:
        ## Calcualte the index size frequencies 
        l,p = indexfreq(indexsizes(loadbam(inpath)),delta=binsize,maxsize=maxindexsize)
        ## Plot the frequency 
        plt.plot(l,p,'-',label=basenobam(inpath),alpha=0.5)
        ## Format into a dataframe and append to our list of rames 
        dfs.append(pd.DataFrame(p,index=l,columns=[basenobam(inpath)]))
    ## Add a legend 
    #plt.legend(bbox_to_anchor=(1.1,1.25))
    ## Modify the x and y ticks
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ## Add an x and y label
    plt.xlabel('Insert Size (bp)',fontsize=12)
    plt.ylabel('Portion of Fragments',fontsize=12)
    ## save the figure if save path was given 
    plt.savefig(outsavepath,dpi=300,bbox_inches='tight') if outsavepath else None 
    ## return the concat dfs
    return pd.concat(dfs,axis=1).reset_index()

## Ftn for parsing chromoeoms 
def getchroms(sam_file,xlist):
    """Ftn for returning a list of chromosomes from a sam file."""
    ## Return the chromosome names 
    return [c for c in sam_file.references if c not in xlist]

## Getting chrom size from sam file
def chromlen(sam_file,chrom):
    """Returns a chromosome length using a bamfile object from pysam given a chromosome name (for e.g. chr1)."""
    ## Return the length given the chromosome name
    return np.array(sam_file.lengths)[(np.array(sam_file.references)==chrom)].min()

## Ftn for binning and counting fragments around chromosomes
def chrombin(sam_file,chrom,window,count=False):
    """Generates windowed regions along a chromosome using a sam file object from pysam.\nThe size of regions is set by WINDOW."""
    ## Gather the sequence length and reset bedcounts
    seqlen, bed_counts = chromlen(sam_file,chrom), []
    ## Reset columns if we need
    colns = ['Chrom','Left','Right','Fragments'] if count else ['Chrom','Left','Right']
    ## Assert the sequence length is greater than zero
    assert (seqlen > 0), "ERROR: The chromosome has a length of zero!"
    ## Set a range given the window size and sequence length
    for j in range(1, seqlen, window):
        ## Set the stop given the window size
        stop = j+window-1 if j+window-1 < seqlen else seqlen
        ## If we are coutning the fragment in this region 
        if count:
            c = len([r for r in sam_file.fetch(chrom,j,stop) if r.is_read1])
            ## Format the counts to append 
            toapp = (chrom, j, stop,c)
        ## Otherwise just append the chrom, left, and right
        else: 
            toapp = (chrom, j, stop)
        ## Append the counts of a given region via the get reads ftn
        bed_counts.append(toapp)   
    ## Format into a dataframe and return 
    return pd.DataFrame(bed_counts,columns=colns)

## Ftn for plotting fragment ranks 
def plotfragmrank(inpaths,resolution=100000,skipchrom=[],pheno = 'Fragments',outsavepath=None,lloc=(0.8,1.25),ax=None):
    """Plots the percent of fragments within genomic bins with respect to the maximum covergae as a function of bin percentage."""
    ## Set the axis if it was give
    plt.sca(ax) if ax else None 
    ## initilze a list of dfs
    dfs = []
    ## Plot a one-to-one line
    plt.plot([0,1],[0,1],'k--',alpha=0.5)
    ## Iterate thru the paths 
    for inpath in inpaths:
        ## Load the bam file 
        bam_file = loadbam(inpath)
        ## Gather the chromlist
        chrlist = getchroms(bam_file,skipchrom)
        ## Calc the fragment counts
        tmp = pd.concat([chrombin(bam_file,c,resolution,count=True) for c in chrlist])
        ## sort the counts
        sort_tmp = tmp[pheno].sort_values().values
        ## Gather the sorted unique values 
        uf = sorted(np.unique(sort_tmp))
        ## Calcualte the percentile 
        newy = np.array([np.max(sort_tmp[(sort_tmp<=u)])/np.max(uf) for u in uf])
        ## Re-calc x
        newx = np.arange(len(uf))/len(uf)
        ## Plot the profile 
        plt.plot(newx,newy,label=basenobam(inpath),alpha=0.5)
        ## Append new df
        dfs.append(pd.DataFrame(newy,index=newx,columns=[basenobam(inpath)]))
    ## Add a legend 
    plt.legend(bbox_to_anchor=lloc)
    ## Modify the x and y ticks
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ## Annotate the x and y label 
    plt.ylabel('Fraction of Fragments\nw.r.t. Highest Coverage',fontsize=12)
    plt.xlabel('Percentile of %s kb Bins'%int(resolution/1000), fontsize=12)
    ## save the figure if save path was given 
    plt.savefig(outsavepath,dpi=300,bbox_inches='tight') if outsavepath else None     
    ## return the concat dfs
    return pd.concat(dfs,axis=1).reset_index()

## ---------------------------------------------- DEFAULT VARIABLES ---------------------------------------------------- ## 
## Set description of this library and scirpt
description = 'Calculates the distribution of fragment lengths from input bam files.'

## Set defaults
windowsize = 100000 
myexcludes = ['chrM']
deltabp = 10 
savepath = './fragment.dist.png'
max_index = 1000 
legendloc = (0.8,1.25)

## Set help messages
b_help = "Path(s) to input BAM files."
W_help = "Size in bp used to construct genomic windows, tiled per chromosome from left to right, for calculating fpkm (default: %s)."%windowsize
X_help = "Names of chromosomes to exclude from analysis (default: %s)."%', '.join(myexcludes)
D_help = "Length (in bp) to merge size of fragments in forming distribution (default: %s). Increasing will smooth distirbution profile."%deltabp
S_help = "Path and name of output diagnostic figure (default: %s)."%savepath
M_help = "The expected maximum size of a fragment (default: %s)."%max_index
L_help = "The location of the legend from the right panel of figure (defualt: (%s,%s))"%legendloc

## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## Load in argparser
    import argparse

    ## Set parser
    parser = argparse.ArgumentParser(description=description)

    ## Add required arguments 
    parser.add_argument("-b", "--bam-files",  dest="b", required=True,  type=str,   help=b_help, nargs='+')
        
    ## Add optional arguments
    parser.add_argument("-W", "--binsize",    dest="W", required=False, type=int,   help=W_help, default=windowsize, metavar='n')
    parser.add_argument("-X", "--exclude",    dest="X", required=False, type=list,  help=X_help, default=myexcludes, metavar='chrM', nargs='+')
    parser.add_argument("-D", "--delta",      dest="D", required=False, type=int,   help=D_help, default=deltabp,    metavar='n')
    parser.add_argument("-S", "--save-path",  dest="S", required=False, type=str,   help=S_help, default=savepath,   metavar='./path/to/png')
    parser.add_argument("-M", "--max-index",  dest="M", required=False, type=int,   help=M_help, default=max_index,  metavar='n')
    parser.add_argument("-L", "--legend-loc", dest="L", required=False, type=tuple, help=L_help, default=legendloc,  metavar='(x,y)')

    ## Parse the arguments
    args = parser.parse_args()

    ## Loadin the correct backend for matplotlib
    import matplotlib

    ## Set the agg backend 
    matplotlib.use('Agg')

    ## Set required vars
    bampaths, windowsize, excludes, deltabp, savepath, maxix, legloc = args.b, args.W, args.X, args.D, args.S, args.M, args.L

    ## Check that we have bam files
    assert len(bampaths) > 0, "ERROR: No bam files were passed to the script."
    
    ## Check all the input bam files are bam files
    for b in bampaths: ## Check if the input bam file is a bam file
        assert isbam(b), "ERROR: The input bam file -- %s -- is not a bamfile!"%b
        assert hasix(b), "ERROR: The input bam file -- %s -- is not indexed!"%b

    ## Call a subplot, set facecolor to white
    fig,ax = plt.subplots(1,2,figsize=(10,4))
    fig.set_facecolor('w')

    ## Plot the fragments
    pfd = plotfragments(bampaths,maxindexsize=maxix,binsize=deltabp,ax=ax[0])
    ## Plot and calcualte binned ranks
    pfr = plotfragmrank(bampaths,resolution=windowsize,skipchrom=excludes,lloc=legloc,ax=ax[1])

    ## Adjust the subplot spacing
    plt.subplots_adjust(wspace=0.25)
    ## Save out the figure 
    plt.savefig(savepath,dpi=300,bbox_inches='tight') if savepath else None 

    ## Save out the dataframes
    if savepath and savepath.split('.')[-1] == 'png':
        ## Set the save path and save out the fragment distribution df 
        pfd_save = savepath.split('.png')[0] + '.frag.dist.csv'
        pfd.to_csv(pfd_save,index=False)
        ## Set the save path and save out the fragment distribution df 
        pfr_save = savepath.split('.png')[0] + '.frag.rank.csv'
        pfr.to_csv(pfr_save,index=False)
## End of file 