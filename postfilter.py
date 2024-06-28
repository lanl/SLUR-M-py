#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
filterdesc = "Filters an input tsv file for properly mapped Hi-C contacts, filtering on mapping quality and intra-chromosomal distance."
## ----------------------------- v 0.0.0 ------------------------------ ## 
##   Post (Hi-C) Filter   
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
## Bring in pandas and np 
import pandas as pd, numpy as np, re 

## Bring in path from slurpy
from defaults import pathexists, juicer_cols, reset, hicsep

## Load in samnames 
from pysamtools import getmatchsum, listzip, writeset, dictzip, loadref, samnames, samdict

## ------------------------ FUNCTION and VARIABLE FILTERING ------------------------------ ## 
## Make list of unique bits
ubits = [2**i for i in range(12)]

## List descripiotn of flags
descr = ['Multiseq',    # template having multiple segments in sequencing
         'Propalign',   # each segment properly aligned according to the aligner
         'Segmunmp',    # segment unmapped
         'Nextsegunmp', # next segment in the template unmapped
         'Seqrev',      # SEQ being reverse complemented
         'Nextseqrev',  # SEQ of the next segment in the template being reverse complemented
         'Read1',       # the first segment in the template
         'Read2',       # the last segment in the template
         'Secalign',    # secondary alignment
         'Nopass',      # not passing filters, such as platform/vendor quality controls
         'Pcrdup',      # PCR or optical duplicate
         'Supalign']    # supplementary alignment
         
## Zip and make into a dict 
bitflags = dict(zip(ubits,descr))

## Check our work
assert len(bitflags) == 12

## Write the cigar strings 
cigarlist = ['M','I','D','N','S','H','P','=','X']

## Format cigar dict
cigardict = dict(zip(cigarlist,range(len(cigarlist))))

## Ftn for spliting cigar
def splitcigar(s,strand):
    c = re.split('(\d+)',s)[1:]
    d = np.array([*map(cigardict.get,c[1::2])],dtype=int)
    n = np.array(c[::2],dtype=int)
    l = listzip(d,n)
    return l if strand else l[::-1]

## Ftn for getting set of of flags from df
def getflags(df):
    ## Return list of flags 
    return list(set(df['Flag']))

## Ftn for adding alignment flags to sam df
def makeflagdf(df):
    ## Gather list of unique flags 
    flags = getflags(df)
    ## Calcualte a dataframe with bit inforamtion for each unqiue flag
    mbits = pd.DataFrame([[f & l for l in ubits] for f in flags],columns=descr)
    ## append the flags
    mbits['Flag'] = list(flags)
    ## Merge the orintations
    #obits = mbits.merge(orientations)
    ## Merge with origional df and return 
    return df.merge(mbits[['Flag','Seqrev','Read1','Read2']]).sort_values('Qname').reset_index(drop=True)

## Ftn for giving back readnames that are miss matched 
def missingmate(df):
    if len(set(df.Qname))*2 != df.shape[0]:
        ## Print a warning 
        print("WARNING: We are missing a mate. Applying secondary filters.")

        ## Find the reads missing a mate (or with more than one)
        ## Take the unique read names and their coutn
        l,c = np.unique(df.Qname,return_counts=True)

        ## Find those with more or less than on mate.
        weirdos = l[(c != 2)]
    else: ## Otherwise return an empty list 
        weirdos = []
    ## Return the set of missing mates 
    return set(weirdos)

## Ftn for appending to df
def appenddf(df,outpath,thesep=hicsep):
    ## Append the df
    df.to_csv(outpath,mode='a' if pathexists(outpath) else 'w',index=False,header= not pathexists(outpath),sep=thesep)
    pass 

## Ftn for formating long dataframe from juicer 
def formatlong(indf,pos_cols0 = ['Seqrev','Rname','Pos','Flag'],pos_cols1 = ['Mapq','Cigar','Seq']):
    """
    Returns a list of mappign values from an input set of rows (representing a Hi-C alignment) to match the juicer format:

    <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <cigar1> <sequence1> <mapq2> <cigar2> <sequence2> <readname1>

    The numeric value of the chormosome is also added.
    """
    ## Set the index as the Qname
    indf.index = indf.Qname
    ## Concat the rows and column spaces into long format 
    long = pd.concat([ indf[pos_cols0][::2],indf[pos_cols0][1::2],
                       indf[pos_cols1][::2],indf[pos_cols1][1::2],
                       indf[['Qname']][::2], 
                       indf[['Chrn']][::2],indf[['Chrn']][1::2],
                       np.abs(indf[['Tlen']][::2])
                       ],axis=1)
    ## Reset column names
    long.columns = juicer_cols + ['Chrn1','Chrn2','Tlen']
    ## Return the longdf
    return long

## Ftn for makign dictionary of chromosome names and id numbers 
def makechromdict(inrefpath):
    ## If we 
    if pathexists(inrefpath+'.fai'):
        ## Load in the chromosome dataframe
        chrdf = pd.read_csv(inrefpath+'.fai',sep='\t',header=None)[[0]]
        ## Format chromosome names into a list 
        chrdict = dictzip(chrdf.T.values[0],chrdf.index)
    else: ## otherwise, load in the ref and make a chromosome dictionary 
        ref = loadref(inrefpath)
        chrdict = dict([(r.name,i) for i,r in enumerate(ref)])
    ## Return the dict
    return chrdict

## Ftn for returning orientation of records 
def getoriented(s1,s2):
    """
    Returns the paired read orientation types, shown below:

    read 1                                  read 2
    -------------->      ...      <---------------
    inward


    <--------------      ...      --------------->
    outward


    <--------------      ...      <---------------
    left 


    -------------->      ...      --------------->
    right

    -----------------------------------------------
    |||||||||||||||||||||||||||||||||||||||||||||||
    --------------------genome---------------------

    """
    ## Gather the strand codes 
    if (s1 == 0) and (s2 ==0):
        orientation = 'Right'
    elif (s1 == 0) and (s2 >0):
        orientation = 'Inward'
    elif (s1 > 0) and (s2 == 0):
        orientation = 'Outward'
    elif (s1 > 0) and (s2 > 0):
        orientation = 'Left'
    else: ## Print an error 
        print("WARNING: We have encountered an orientation we did not recognize")
        orientation = 'Unknown'
    return orientation

## Ftn for post filtering recroeds 
def postfilter(inmapping,outdfpath,setdistance,mapq,chrdict):
    ## add mapped flags
    mmapped = makeflagdf(inmapping)
    ## Check our work
    assert mmapped.shape[0] == inmapping.shape[0], "ERROR: There was a critical error in handeling read pairs."
    
    ## Gather the unique cigars from positivly orinted reads, calcualte and add the matcsum 
    poscigars = mmapped[(mmapped.Seqrev==0)][['Cigar']].drop_duplicates().copy().reset_index(drop=True)
    poscigars['Matchsum'] = poscigars.Cigar.apply(splitcigar,args=[1]).apply(getmatchsum)

    ## Gather the unique cigars form negativel orinted reads, calcualte and add the matcsum 
    negcigars = mmapped[(mmapped.Seqrev>0)][['Cigar']].drop_duplicates().copy().reset_index(drop=True)
    negcigars['Matchsum'] = negcigars.Cigar.apply(splitcigar,args=[0]).apply(getmatchsum)

    ## Calculate the number of positively and negativel oriented reads 
    pos_shape = mmapped[(mmapped.Seqrev==0)].shape[0]
    neg_shape = mmapped[(mmapped.Seqrev>0) ].shape[0]
    
    ## Merge the cigar values 
    pos = mmapped[(mmapped.Seqrev==0)].merge(poscigars)
    neg = mmapped[(mmapped.Seqrev>0) ].merge(negcigars)

    ## Check our shapes 
    assert (pos.shape[0] == pos_shape) and (neg.shape[0] == neg_shape)

    ## make a new map with the cigar values 
    newmap = pd.concat([pos,neg]).sort_values('Qname').reset_index(drop=True)
    ## Make the matsum column an int
    newmap['Matchsum'] = newmap['Matchsum'].apply(int)

    ## Check our shape and readnames
    assert newmap.shape[0] == mmapped.shape[0]
    assert not len(set(newmap.Qname) - set(inmapping.Qname))
    
    ## Gather the minimized index for the match sum column and set to list 
    r1_ix = newmap[(newmap.Read1>0)].groupby('Qname').Matchsum.idxmin().tolist()
    r2_ix = newmap[(newmap.Read2>0)].groupby('Qname').Matchsum.idxmin().tolist()

    ## Gather the pairs for each read
    pairmap = newmap.loc[sorted(r1_ix+r2_ix),:]
    
    ## Check the pairing
    oddreads = missingmate(pairmap)
    
    ## Quality filters, gather the set of failigin reads
    lowqual_reads = set(pairmap[(pairmap.Mapq<mapq)].Qname)
    
    ## apply distance threshodling on intra-chromosomal contacts
    #small_tlen = set(pairmap[(pairmap.Rnext=='=')&(pairmap.Tlen.abs()<=setdistance)].Qname) if setdistance else set()
    
    ## Filter the pair map 
    #pairmap = pairmap[~(pairmap.Qname.isin(small_tlen | lowqual_reads | oddreads))] 
    pairmap = pairmap[~(pairmap.Qname.isin(lowqual_reads | oddreads))] 
    ## Add the chromosome as integer
    pairmap['Chrn'] = pairmap.Rname.replace(chrdict)
    ## Sort the mapings by the Qname, chrn and 
    pairmap = pairmap.sort_values(['Qname','Chrn','Pos']).reset_index(drop=True)

    ## Make the long df
    longtmp = formatlong(pairmap)
    longtmp['Tlen'] = (longtmp.Pos1 - longtmp.Pos2).abs()

    ## Gather the read pairs with smaller mapped distance than threshold 
    small_tlen = set(longtmp[(longtmp.Tlen<=setdistance)].Qname1)

    ## Format the long 
    long = longtmp[~(longtmp.Qname1.isin(small_tlen))].copy()
    
    ## Add the orientation
    long["Orientation"] = long.apply(lambda x: getoriented(x["Str1"], x["Str2"]), axis=1) if (long.shape[0] > 0) else None

    ## Save out the pairs and pass back to main
    appenddf(long,outdfpath) if long.shape[0] else None

    ## Return the set of removed read names
    return oddreads, lowqual_reads, small_tlen

## Set help strings
T_help = "Path to input .sam file or tsv file from bwa."
S_help = "Extension to seperate from input file; used to generate output file names (default: .txt)."

## Load in defaults from slurpy 
from defaults import chunks, map_q_thres, Q_help, Z_help, set_distance, D_help, r_help 

## ---------------------------- MAIN SCRIPT & ARGUMENT PARSING -------------------------- ## 
## If the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = filterdesc)

    ## Add the required arguments
    parser.add_argument("-t", dest="T", type=str,  required=True,  help=T_help, metavar='./path/to/input.sam'        )
    parser.add_argument("-r", dest="R", type=str,  required=True,  help=r_help, metavar='./path/to/ref.fai'          )
    
    ## Add optional variables 
    parser.add_argument("-Q", dest="Q", type=int,  required=False, help=Q_help, metavar='n',    default=map_q_thres  )
    parser.add_argument("-D", dest="D", type=int,  required=False, help=D_help, metavar='n',    default=set_distance )
    parser.add_argument("-Z", dest="Z", type=int,  required=False, help=Z_help, metavar='n',    default=chunks       )
    parser.add_argument("-S", dest="S", type=str,  required=False, help=S_help, metavar='.txt', default='.txt'       )

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set inputs
    inpath       = inputs.T  ## Path to input same file 
    refpath      = inputs.R  ## input path to reference 
    map_q        = inputs.Q  ## map quality threshold 
    set_distance = inputs.D  ## Set linear distance for filtering 
    chunk_size   = inputs.Z  ## Chunk-size
    exten        = inputs.S  ## the file extension used 

    ## Load in file
    filehead = inpath.split(exten)[0]

    ## Set file ends 
    fileends = ['oddling','lowqual','distance','filtered']

    ## Set file ends and new fiel namees
    newfiles = [filehead+'.%s.txt'%f for f in fileends] 

    ## Reset, removing previous runs
    reset(newfiles)
    
    ## Make into a dict 
    savefiledict = dictzip(fileends,newfiles)

    ## Format the chromosome dict
    chromdict = makechromdict(refpath)

    ## Initilzse lastc and sets of failed reads
    lastc = []
    odd_reads = set()
    low_quals = set()
    distancef = set()

    ## Open the inpath 
    with pd.read_csv(inpath,chunksize=chunk_size,sep='\t',usecols=samnames,dtype=samdict) as samchunks:
        ## Iteraate thru the chunks 
        for inmap in samchunks:
            ## Print the shape

            ## If this is the start of our second time or more thru this loop
            inmap = pd.concat([lastc,inmap],axis=0) if len(lastc) else inmap
                
            ## Gather the last Q name 
            lastq = inmap.tail(3).Qname.tolist()[-1]
            
            ## parse off the last Q
            lastc = inmap[(inmap.Qname==lastq)].copy()
            
            ## Cut off the end if this chunk is the smallest and their for the last
            tomap = inmap[(inmap.Qname!=lastq)]

            ## Call our filtering ftn
            oddlings, belowqual, fdistance = postfilter(tomap,savefiledict['filtered'],set_distance,map_q,chromdict)

            ## Append the set or failed read naems
            odd_reads = odd_reads | oddlings
            low_quals = low_quals | belowqual
            distancef = distancef | fdistance
            
        ## Make sure we get the last bit if it wasn't included, run the last chunk. I don't like this coding and will retrun to it later. 
        oddlings, belowqual, fdistance = postfilter(lastc,savefiledict['filtered'],set_distance,map_q,chromdict) if (not (lastq in set(tomap.Qname)) and len(lastc)) else (set(), set(), set())

        ## Append the set or failed read naems
        odd_reads = odd_reads | oddlings
        low_quals = low_quals | belowqual
        distancef = distancef | fdistance

        ## Write the failed reads out to file
        [(writeset(savefiledict[k],theset) if len(theset) else None) for k,theset in zip(['oddling','lowqual','distance'],[odd_reads,low_quals,distancef])]
## End of file 