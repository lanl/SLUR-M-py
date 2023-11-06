#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
####################################################################
##    COUNTBAMS - Counts reads within bam files made by SLURPY    ##
####################################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ----------------------------------------- LOAD IN MODULES -------------------------------------------- ## 
## Load in system, json, and pandas 
import sys, json, pandas as pd 

## Load in check sam 
from pysamtools import loadbam, readset

## Load in sort-glob 
from slurpy import sortglob, fileexists, aligndir, diagdir, debugdir, basenobam

## ---------------------------------------- Variable Setting -------------------------------------------- ##
## Set the count columns, colors vector, and marked ending
countcolumns, mycolors = ['Filename','Mapping','Readcounts'], ['red','blue','pink','green','orange','black','purple','grey']

## -------------------------------------- Function Defining -------------------------------------------- ##
## Ftn for testing if an input string is equal to marked
def ismarked(m):
    """Returns the boolean that m is marked."""
    ## Return negative bool 
    return m == 'marked'

## Write ftn for splitting file name
#def splitfilename(inbam):
#    """Splits an input bam file name."""
#    ## Return the split file 
#    return inbam.split('/')[-1].split('.bam')[0]

## Write ftn for gathering bam map names
def bamnames(inbam):
    """Parses the endding bam file name."""
    ## Return the parsed . joined name 
    return basenobam(inbam).split('.')[-1]

## Write ftn for calculating the number of records
def nrecords(inbam):
    """Calculates the number of records within an input bamfile."""
    ## Return the count
    return len(readset(loadbam(inbam)))

## Ftn to count duplicates
def countdups(inbampath):
    """Counts the unique read names if they are duplicates."""
    ## Calcualte the length of the set 
    return len(set([r.qname for r in loadbam(inbampath) if r.is_duplicate]))

## Write ftn for counting lines in file
def countlines(inpath):
    """Opens a file for reading only and counts the lines in file."""
    ## Set the count variable 
    count = 0 
    ## If we can find the input path
    if fileexists(inpath):
        ## with open, read in the lines 
        with open(inpath,'r') as openfile:
            for line in openfile:
                ## For each line add to the count 
                count += 1
            ## Close the file 
            openfile.close()
    else: ## Otherwise do nothing 
        pass 
    ## Return count
    return count 

## Ftn for returning unmapped reads
def getunmapped(inpath):
    """Returns counts of the unmapped reads."""
    ## Return the split informations 
    return inpath.split('/')[-1].split('.txt')[0],'unmapped',countlines(inpath)

## Ftn for iteratively counting bam 
def calcbamcounts(inbampaths,unmappedpaths):
    """Iterates thru a list of input bam paths and calculates the number of unique read names."""
    ## Initilize bamcounts 
    bamcounts = []
    ## Iterate thru the input bams
    for inbampath in inbampaths:
        ## Calcualte the map name 
        mapname = bamnames(inbampath)
        ## if this is the marked bam file name, we count the duplicates, otherwise append the counts, append the results 
        bamcounts.append((basenobam(inbampath),'duplicates' if ismarked(mapname) else mapname,countdups(inbampath) if ismarked(mapname) else nrecords(inbampath)))
    ## append the unmapped counts
    for unmpath in unmappedpaths:
        bamcounts.append(getunmapped(unmpath))
    ## Return the bamcounts as a dataframe 
    return pd.DataFrame(bamcounts,columns = countcolumns)

## Ftn for calculating txt counts
def calctxtcounts(intxts):
    ## Initilize bamcounts 
    bamcounts = []
    ## Iterate thru
    for inpath in intxts:
        ## If the extension is a txt file 
        if inpath.split('.')[-1] == 'txt':
            ## Gather the second to last by '.'
            mapname = inpath.split('.txt')[0].split('.')[-1]
        ## If it is a bedpe file 
        elif inpath.split('.')[-1] == 'bedpe':
            ## make it a hic 
            mapname = 'hic'
        else: ## Otherwise take the end 
            mapname = inpath.split('.')[-1]
        ## Iterate thru 
        bamcounts.append((basenobam(inpath),mapname,countlines(inpath)))
    ## Return the dataframe
    return pd.DataFrame(bamcounts,columns = countcolumns)

## Ftn for parsing json file
def parsejson(injson):
    """Load and parses an input json file, returns object as json dictionary."""
    ## Loads the json file
    with open(injson,'r') as infile:
        data = json.load(infile)
    ## Return the parsed data
    return data

## Ftn for getting read totals 
def jsontotals(injson,beforeorafter):
    """Parses a json dictionary from fastp to gather read counts. Beforeorafter is a boolean controlling return of counts befor (TRUE) or after filtering (FALSE)."""
    ## Return the parsed jason
    return int(parsejson(injson)['summary']['before_filtering' if beforeorafter else 'after_filtering']['total_reads'])

## ------------------------------------------------------ BODY of MAIN EXECUTABLE --------------------------------------------------------- ##
## If the script is called 
if __name__ == "__main__":
    ## Gather the run name from system input 
    run_name, expmode = sys.argv[1], sys.argv[2]

    ## if we are in hicmode
    if expmode == 'hic':
        ## Gather the bampe file, and the the txt files, including the unmapped
        hic_files = sortglob(f'./{aligndir}/*.bedpe') + sortglob(f'./{aligndir}/*.txt')
        ## Check whtat we have files
        assert len(hic_files), "ERROR: Unable to detect hic files."

        ## Initilzse the bam counts
        thebamcounts = calctxtcounts(hic_files)
        
    else: ## Set the bam paths, initilize counts
        bampaths = sortglob(f'./{aligndir}/*.bam')
        ## Check that we have bam files
        assert len(bampaths), "ERROR: Unable to detect bam files on this given wild card path: ./aligned/*.bam"

        ## count the unmapped txt, gather the counts
        unmapped_paths = sortglob(f'./{aligndir}/*.unmapped.txt')
        ## Check our work
        assert len(unmapped_paths), "ERROR: Unable to detect list of unmapped read names!"

        ## Iteratively calculate the n-records
        thebamcounts = calcbamcounts(bampaths,unmapped_paths)

    ## Set the sample name for each count
    thebamcounts['Sample'] = [s.split('.')[0] for s in thebamcounts.Filename.tolist()]

    ## Gather the unique smaples
    sample_names = thebamcounts.Sample.unique()

    ## Set the json itertion 
    i = 1 if expmode == 'hic' else 0

    ## Load in json files, set the json paths 
    try_json_paths = [f'./{debugdir}/{s}.fastp.log.txt.{i}.json' for s in sample_names]

    ## Iterate thru the json paths
    for samplename,json_path in zip(sample_names,try_json_paths):
        ## If the file exists 
        if fileexists(json_path):

            ## Gather the starting total and count after filtering 
            true_total, after_filt = jsontotals(json_path,True), jsontotals(json_path,False)

            ## Sum the read coutns for this sample, and double them to pairs
            sample_map = 2*thebamcounts[(thebamcounts.Sample==samplename) & ~(thebamcounts.Mapping.isin(['total','q30','duplicates']))].Readcounts.sum()
            ## Check the numbers of counts match
            if not (sample_map == after_filt):
                print(f'WARNING: The counts of reads post fastp filtering ({after_filt}) and post mapping ({sample_map}) do not match!')

            ## Calc the low quality reads 
            thebamcounts.loc[len(thebamcounts.index)] = ['', 'lowqual', (true_total-after_filt)/2,samplename] 
            ## Append the new entry
            thebamcounts.loc[len(thebamcounts.index)] = ['', 'total', true_total/2,samplename] 
            
        ## Otherwise print a wrning 
        else:
            print('WARNING: Unable to find json file: %s'%json_path)

    ## Format the bamcounts into a pandas df and write to file 
    thebamcounts.to_csv(f'./{diagdir}/{run_name}.bam.counts.csv',index=False)

    ## Convert the counts to million pairs
    thebamcounts['Paired Reads ( Millions )'] = thebamcounts.Readcounts/(10**6)

    ## Make a color dictionary 
    colordict = dict(zip(thebamcounts.Mapping.unique(),mycolors))

    ## Load in matplot lib 
    from matplotlib import pyplot as plt

    ## Bring in seaborn
    import seaborn as sns 

    ## Call a figure, set facecolor
    fig, ax = plt.subplots(1,1,figsize=(7,6))
    fig.set_facecolor('w')

    ## Plot the barplot 
    sns.barplot(y='Sample',x='Paired Reads ( Millions )',hue='Mapping',data=thebamcounts,palette=colordict,ax=ax)
    
    ## Saveout the figure
    plt.savefig(f'./{diagdir}/{run_name}.bam.counts.png',dpi=300,bbox_inches='tight')
## End of file