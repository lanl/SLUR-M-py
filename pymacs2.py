#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
########################################
###    FROMATTING MACS2 FUNCTIONS    ###
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
## ------------------------------------------------------------- MACS2 Functions ----------------------------------------------------------------- ## 
## Load in pandas 
import pandas as pd 

## Load in scripts dir 
from defaults import scriptsdir

## Ftn for formating length parameter
def formatlen(minlen):
    """Formats and returns the max length flag and parameter for a call to macs2."""
    ## Return the len
    return '--min-length %s'%minlen if minlen else ''

## Ftn for formating gap parameter
def formatgap(maxgap):
    """Formats and returns the max gap flag and parameter for a call to macs2."""
    ## Return the gap
    return '--max-gap %s'%maxgap if maxgap else ''

## Ftn for making input into a list of inputs
def makelist(input):
    """Makes an input variable into a list."""
    ## Return the listed input
    return input if (type(input) == list) else [input]

## Ftn for formating inputs
def sjoin(input):
    """Joins an input into a list, seperated by a space."""
    ## Make the input a list 
    inputs = makelist(input)
    ## Make the inputs into a joined list 
    return inputs[0] if (len(inputs) == 1) else ' '.join(inputs)

## Ftn for formatting control
def formatcontrol(incontrols):
    """Formats the input controls for a call to macs2 in chip mode."""
    ## Return the formated input control 
    return f'-c {sjoin(incontrols)}' if incontrols else ''

## Ftn for formatting inputs
def formatinput(inputs):
    """Formats the input bam files for a call to macs2."""    
    ## Return the joined inputs for macs2
    return f'-t {sjoin(inputs)}' if inputs else ''

## Set options for atac-seq opts
atacopts = '--keep-dup all --call-summits -B --SPMR --nolambda'

## Set options for chip-seq opts
narrow_chip = '--keep-dup all --call-summits -B --SPMR'
broad_chip  = '--keep-dup all --broad -B --SPMR'

## Write ftn for calling macs2 with atac seq data
def peakattack(inbams,n,report,broad=False,outdir='./macs2',gsize='hs',mg=None,ml=None,extraoptions=None,incontrols=None):
    """Envokes macs2 callpeak function on for a run of the slurpy pipeline (n) on input bam file in bampe mode using the input genome size (g), maximum gap (ml), and minimum peak length (ml)."""
    ## Format the input options, if the options are set explicitly, and if the input control is set for chip experiment 
    if incontrols and broad:
        opts = broad_chip
    elif incontrols and not broad:
        opts = narrow_chip
    else: ## otherwise this is an atac exp 
        opts = atacopts
    ## Add the additional optsions 
    opts = opts + (' ' + extraoptions if extraoptions else '')
    ## Format the macs2 callpeak command
    return [f'macs2 callpeak {formatinput(inbams)} {formatcontrol(incontrols)} -n {n} -g {gsize} -f BAMPE --outdir {outdir} {opts} {formatgap(mg)} {formatlen(ml)} 2>> {report}\n', f'{scriptsdir}/myecho.py Finished calling peaks in {sjoin(inbams)} with macs2 >> {report}\n']

## Set the narrow peak names
peaknames = ['Chrom','Start','End','Name','Score','Strand','Fold_change','-log10pvalue','-log10qvalue','Sumpos']

## Write ftn for loading in narrow peak file
def loadnarrowpeak(path,peakcols = peaknames,sep='\t'):
    """Loads in a narrow peak file from MACS2."""
    ## Load in narrow peak file and return 
    return pd.read_csv(path,sep=sep,names=peakcols)    

## Wrt ftn for slimming down narrow peak file
def slimdown(inbed,coi=['Chrom','Start','End']):
    """Returns a slimmer bed file from input, returning only columns of interst (COI) and droping duplicate rows."""
    ## Return only the columns of interest
    return inbed[coi].drop_duplicates()

## ---------------------------------------- Write ftns for peak analysis ------------------------------------------------ ## 
## Ftn for calc reads in peaks 
def readsinpeaks(inbam,peaks):
    """Summs the counts of pair-reads within an input peak file."""
    ## Return the counts 
    return sum([len([ r for r in inbam.fetch(contig=row.Chrom,start=row.Start,stop=row.End) if r.is_read1 and not r.is_supplementary]) for i,row in peaks.iterrows()])

## Ftn for calculating unique reads in a bam file 
def readsinbam(inbam):
    """Gathers read names from an input bam file and calculates the unique set of read names."""
    ## Return the length of the set of read names 
    return len(set([r.qname for r in inbam]))

## Ftn for calculating bam map from chromosome 
def makemap(inbam):
    """Generates a chromosome map from an input bam file."""
    ## Return the dataframe 
    return pd.DataFrame(zip(inbam.references,inbam.lengths),columns=['Chrom','Length'])

## ---------------------------------------------- DEFAULT VARIABLES ---------------------------------------------------- ## 
## Set description of this library and scirpt
description = 'Calculates the fraction of reads within peaks from input bam and bed files.'

## Set defaults
savename, dplace = './frip.stats.csv', 4

## Set help messages
b_help = "Path(s) to input BAM files."
p_help = "Path(s) to input peak (BED) files from macs2."
s_help = "Path and name of output diagnostic statistics."
d_help = "Decimal place used to calcualte and save statistics (Default: %s)."%dplace

## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## ------------------------------------------- MODULE LOADING ---------------------------------------------------- ## 
    ## Load in pandas and arg parser
    import argparse 

    ## Ftn for parsing files 
    from defaults import sortglob, aligndir, macs2dir, diagdir

    ## Load bam ftn 
    from pysamtools import loadbam, isbam, hasix

    ## ------------------------------------------ PARSER SETTING ---------------------------------------------------- ## 
    ## Set parser
    parser = argparse.ArgumentParser(description=description)

    ## Add optional arguments
    parser.add_argument("-b", "--bam-files",  dest="b", required=False, type=list, help=b_help, default=[],nargs='+')
    parser.add_argument("-p", "--bed-files",  dest="p", required=False, type=list, help=p_help, default=[],nargs='+')
    parser.add_argument("-s", "--save-path",  dest="s", required=False, type=str,  help=s_help, default=None        )  
    parser.add_argument("-d", "--decimals",   dest="d", required=False, type=int,  help=d_help, default=dplace      )

    ## Parse the arguments
    args = parser.parse_args()

    ## ---------------------------------------- VARIABLE SETTING ---------------------------------------------------- ## 
    ## Gather inputs 
    bams_paths, peak_paths, save_path, dplace = args.b, args.p, args.s, args.d

    ## Re-assign bam paths form the local aligned dir by gathering the primary mapped bam files from the aligned dir if run name was passed 
    bams_paths = bams_paths if len(bams_paths) else sortglob(f'./{aligndir}/*.primary.*.bam') 

    ## Re-assign bed paths from local macs2 dir if we were passed a run name, gather the peaks from the macs2 dir
    peak_paths = peak_paths if len(peak_paths) else sortglob(f'./{macs2dir}/*_peaks.*Peak') 

    ## Reset the save path if none was given 
    save_path = save_path if save_path else f'./{diagdir}/{savename}' 

    ## Check we have paths 
    assert len(bams_paths), "ERROR: No bam files were detected!"
    assert len(peak_paths), "ERROR: No macs2 peak files were detected!"

    ## Check all the input bam files are bam files
    for b in bams_paths: ## Check if the input bam file is a bam file
        assert isbam(b), "ERROR: The input bam file -- %s -- is not a bamfile!"%b
        assert hasix(b), "ERROR: The input bam file -- %s -- is not indexed!"%b

    ## --------------------------------------- FRiP Calculation --------------------------------------------------- ## 
    ## Initlizse peak count and info
    peak_info = []
    ## Iterate thru the narrow peak files 
    for peak_path in peak_paths:
        ## Load in the narrow peak file and initilize the read counts and totals
        narrow, read_total, read_peaks  = loadnarrowpeak(peak_path), [], []
        ## Slim down the narrow peak file
        slimpeak = slimdown(narrow)
        ## Calculate the total bp covered
        total_bp = sum(slimpeak['End'].values - slimpeak['Start'].values)
        
        ## Iterate thru the bam files 
        for bam_path in bams_paths:
            ## Load in the bam
            bam = loadbam(bam_path)
            ## Count reads in bam file
            read_total.append(readsinbam(bam))
            ## Count reads in peak
            read_peaks.append(readsinpeaks(bam,slimpeak))
            
        ## Calc the frip 
        peak_frip = round(sum(read_peaks)/sum(read_total),dplace)
        ## Append the resutls
        peak_info.append([peak_path.split('/')[-1],narrow.shape[0],slimpeak.shape[0],peak_frip,total_bp])

    ## Get the chrom map
    chrommap = makemap(bam)
    ## Gather the chromosome list
    chrlist = slimpeak.Chrom.unique()
    ## Calcualte the genome size 
    genomesize = chrommap[(chrommap.Chrom.isin(chrlist))].Length.sum()
        
    ## Format into a df
    peak_info = pd.DataFrame(peak_info,columns = ['Peak File','Summits','Peaks','FRiP','BP'])
    ## Calculate the perecnt genome
    peak_info['Percent'] = round(100*peak_info.BP/genomesize,dplace)
    ## Save the peak info
    peak_info.to_csv(save_path,index=False)
## End of file 