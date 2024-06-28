#!/usr/bin/env python
debug = False 
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
########################################
###    Hi-C Compartment Analysis     ###
########################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ------------- Function Defining ------------------ ##
## Import fanc, sys, nupy and and padns 
import fanc, numpy as np, pandas as pd

## Load in narrow peak file 
from pymacs2 import loadnarrowpeak

## Bring in matplot lib 
from matplotlib import pyplot as plt

## Load in statmodels 
from statsmodels.api import nonparametric

## Returns df by chormosmoe 
def bychrom(df:pd.core.frame.DataFrame, coi:str) -> pd.core.frame.DataFrame:
    ## Return the chrom data 
    return df[(df.Chrom==coi)]

## Ftn for getting bin coutns of peaks 
def bincount(bins:list, df:pd.core.frame.DataFrame) -> np.array:
    ## Counts peaks per bins in peak dataframe 
    return np.array([df[(df.Start>=b) & (df.End<=bins[i+1])].shape[0] for i,b in enumerate(bins[:-1])])

## Ftn for chetting chormosome eigen vector 
def getchrev(evector:list,cdict:dict,coi:str) -> np.array:
    ## Return the eigen vector 
    return np.array(evector[cdict[coi][0]: cdict[coi][1]])

## Ftn for getting sign of array 
def getsign(n: float) -> int:
    ## Return the integer sign 
    return 1 if (n>=0) else -1

## Ftn for turnning off axis spins
def spineoff(az):
    ## Turnoff the spines 
    [az.spines[t].set_visible(False) for t in ('top','right')]
    pass 

## Set column names 
chrom_cols = ['Chrom','Start','End','Eigenvalue','Peakcounts']
    
## Set the description
description = 'Estimates A and B compartmetns from Hi-C and peak called dataframe.'

## Set defaults
fold_chan = 2
bins_size = 250000
ncol = 7
myfs = 12
mydpi = 300
myfigsize = (15,15)
span_size = 15
lowess_fraction = 0.4

## Set help messages
i_help = 'Path to input .hic file.'
p_help = 'Path to input narrow peak file from macs2.'
b_help = 'The binsize or resolution (bp) to use in analysis (default: %s).'%bins_size
f_help = 'The minimum fold change to filter peaks (default: %s).'%fold_chan
o_help = 'Output direcotory to save resluts in csv file to.'
s_help = 'Size of vector to take from front end of lowess regression for calculating correction (default: %s).'%span_size
l_help = 'Between 0 and 1. The fraction of the data used when estimating each y-value (default: %s).'%lowess_fraction  
P_help = 'Boolean flag to skip plotting compartment scores.'

## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## ------------------------------------------- MODULE LOADING ---------------------------------------------------- ## 
    ## Load in pandas and arg parser
    import argparse 
    ## ------------------------------------------ PARSER SETTING ---------------------------------------------------- ## 
    ## Set parser
    parser = argparse.ArgumentParser(description=description)

    ## Set input arguments
    parser.add_argument("-i", "--inhic-file",        dest="i", required=True,  type=str,   help=i_help, metavar= './path/to/in.hic') 
    parser.add_argument("-p", "--peaks-file",        dest="p", required=True,  type=str,   help=p_help, metavar= './path/to/narrow.peaks')

    ## Set default arguments 
    parser.add_argument("-b", "--binsize",           dest="b", required=False, type=int,   help=b_help,  default= bins_size)
    parser.add_argument("-f", "--fold-change",       dest="f", required=False, type=int,   help=f_help,  default= fold_chan)
    parser.add_argument("-o", "--outputpath",        dest="o", required=False, type=str,   help=o_help,  default= None,           metavar = './path/to/saveout.csv')
    parser.add_argument("-s", "--span-size",         dest="s", required=False, type=int,   help=s_help,  default=span_size,       metavar='n')
    parser.add_argument("-l", "--lowess-fraction",   dest="l", required=False, type=float, help=l_help,  default=lowess_fraction, metavar='p')

    ## Set boolean vars
    parser.add_argument("--no-plot",                 dest="P",  help = P_help, action = 'store_true')

    ## Parse the arguments
    args = parser.parse_args()

    ## Set input arguments 
    hicf_path = args.i         ## Path to hic path
    peak_path = args.p         ## Path to narrow peak file 
    bins_size = args.b         ## Set the bin size used in analysis 
    fold_chan = args.f         ## Calc the fold change 
    head_path = args.o         ## gather the head path 
    toplot    = not args.P     ## Set plotting boolean 
    span_s    = args.s         ## Set span of lowess
    lf        = args.l         ## Lowess fraction    

    ## Set the head pathif none was given
    head_path = head_path if head_path else ('/'.join(hicf_path.split('/')[:-1]) + '/')

    ## Reset the head paths last character 
    head_path = head_path if (head_path[-1] == '/') else (head_path + '/')

    ## Check that this path is relative
    head_path = head_path if (head_path[0] == '.') else ('.'+head_path)

    ## Set the savepath
    savepath = head_path + hicf_path.split('/')[-1].split('.hic')[0] + '.compartment.scores.%s.csv'%bins_size
    saveppng = savepath.split('.csv')[0] + '.pdf'

    ## Print the paths 
    print(savepath) if debug else None
    print(saveppng) if debug else None 

    ## Load in peaks
    peaks = loadnarrowpeak(peak_path)

    ## Set the bin str
    binstr = '@%skb'%int(bins_size/1000)

    ## Load in hic file 
    hic_1mb = fanc.load(hicf_path+binstr)

    ## Gather chromosome bins 
    chrlist = hic_1mb.chromosomes()

    ## Set the bin dict of chromosomes
    chrbins = hic_1mb.chromosome_bins

    ## Calcualte a and b compartment 
    ab = fanc.ABCompartmentMatrix.from_hic(hic_1mb)

    ## Call ev ftn for gather eigen vectors
    ev = ab.eigenvector()

    ## Set chrom list
    chrom_dfs = []

    ## If we are plotting 
    if toplot: ## Call a figure, set facecolor 
        ## Set rows
        row_counts = int(np.ceil(len(chrlist)/ncol))
        ## Call a figure
        fig,ax = plt.subplots(row_counts,ncol,figsize=myfigsize,sharex=True,sharey='row')
        fig.set_facecolor('w')
        ## Make a vecotr of rows
        axis_ixs = list(range(int(row_counts*ncol)))

    ## Iterate thru the chromosomes
    for cix,chrom in enumerate(chrlist): 
        ## Filter the peaks for this chromosome 
        filt = bychrom(peaks[(peaks.Fold_change>fold_chan)],chrom)[['Start','End']].drop_duplicates()

        ## Gather the bin bounds
        a,b = hic_1mb.chromosome_bins[chrom]

        ## Gather bins 
        bins = [0] + list(np.cumsum(np.repeat(bins_size,b - a)))

        ## Calc bin counts from df
        bin_counts = bincount(bins,filt)

        ## Calc chrom ev
        chrom_ev = getchrev(ev,chrbins,chrom)

        ## Calculate lengths 
        bc_l = len(bin_counts)
        ce_l = len(chrom_ev)

        ## Check our work, and that the bin sizes match 
        assert (bc_l == ce_l), "ERROR: Bin sizes (%s,%s) do not match for chromosome: %s!"%(bc_l,ce_l,chrom)

        ## Compute a lowess smoothing of the data
        smoothed = nonparametric.lowess(exog=bin_counts, endog=chrom_ev, frac=lf)

        ## Calcualte the correction span
        correction_span = smoothed[:,1][:int(span_s)].mean()

        ## Get the sign 
        thesign = getsign(correction_span)

        ## Set correction method, Do we need to correct this?
        to_correct, iscor = (-1,'\n( Corrected )') if (thesign > 0) else (1,'')

        ## Correct the estimates of open or closed
        corrected = chrom_ev * to_correct

        ## Set the chromosome df
        chrom_df = pd.DataFrame([np.repeat(chrom,len(corrected)),bins[:-1],bins[1:],corrected,bin_counts],index=chrom_cols).T

        ## Append to list 
        chrom_dfs.append(chrom_df)

        ## If we are plotting 
        if toplot:
            ## Set the axis handl
            az = ax.ravel()[cix]
            ## Set the fig axis for this subplot 
            plt.sca(az)

            ## Plot the eigen values 
            plt.plot(bin_counts[(corrected>=0)], corrected[(corrected>=0)], '.', color='tab:red' , rasterized=True)
            plt.plot(bin_counts[(corrected <0)], corrected[(corrected <0)], '.', color='tab:blue', rasterized=True)

            ## Plot the smoothed
            plt.plot(smoothed[:,0], smoothed[:,1]*to_correct, color='k', rasterized=True)
            
            ## Add title to subplot 
            plt.title(chrom + iscor,fontsize=myfs,y=0.8)
            ## Turn off the spines
            spineoff(az)

    ## Remove the last axis 
    if toplot:
        for i in axis_ixs[len(chrlist):]:
            plt.sca(ax.ravel()[i])
            plt.axis('off')

    ## Merge the chrom dfs into genom
    genomic_df = pd.concat(chrom_dfs)

    ## Save-out the dataframe
    genomic_df.to_csv(savepath,index=False)

    ## Annotate the x and yaxis of the plot
    if toplot:
        ## Label x and y 
        fig.text(x=0.5,y=0.07,s='Number of Peaks ' + binstr,fontsize=myfs,va='center',ha='center')
        fig.text(x=0.05,y=0.5,s='Eigen Value',fontsize=myfs,va='center',ha='center',rotation=90)

    ## Save the figure
    plt.savefig(saveppng,dpi=mydpi,bbox_inches='tight') if toplot else None 
## End of file 