#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## ----------------------------- v 8.0.0 ------------------------------ ## 
##   PYSAMTOOLS: Python functions for Working with SAMTOOLS (v 1.5.1)   ##
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
## Bring in subprocess, pysam 
import subprocess, pysam 

## Load in is file ftn 
from os.path import isfile

## ----------------------------------- GENERAL FUNCTIONS --------------------------------- ##
## Turns an input into a list 
def makelist(input):
    """Transforms an input into a list."""
    ## Returns list 
    return input if (type(input) is list) else [input]

## Ftn for making a dicitonary for zipped lists
def dictzip(a,b):
    """Makes a dictionary from zipped items in lists a and b."""
    ## Return the zipped dicts
    return dict(zip(a,b))
    
## Write a ftn for printing given a condition
def ifprint(message,bool):
    """Prints message given the boolean state."""
    ## Print the message
    print(message) if bool else None
    ## Return the boolean
    return bool 

## Ftn for comparing the lengths of two arrays
def lencomp(a,b):
    """Checks if the lengths of lists a and b are equal."""
    ## Check that a and b have the same length
    return (len(a) - len(b)) == 0

## Cuts on the bam extension
def splitbam(inbam):
    """Splits an input text on the string .bam."""
    ## Return the split on the sting
    return inbam.split('.bam')[0]

## Set ftn for making bam file names
def outnames(inbam,mito):
    """Format and return the name of output txt, bam, and bedpe files."""
    ## Format the name of output files 
    name_out_bam     = splitbam(inbam) + '.mapped.bam'    ## The output bam name
    name_placed_txt  = splitbam(inbam) + '.placed.txt'    ## The placed reads (one mate mapped) bam file  
    name_mito_txt    = splitbam(inbam) + '.%s.txt'%mito   ## The mitochondrial mapped reads
    name_unmap_txt   = splitbam(inbam) + '.unmapped.txt'  ## The unmapped read pairs txt file
    name_out_bedpe   = splitbam(inbam) + '.bedpe'         ## The bedpe file 
    ## Return names
    return name_out_bam, name_placed_txt, name_mito_txt, name_unmap_txt, name_out_bedpe

## Ftn for returning split sub names
def splitsubnames(mito):
    """Returns a formated list of middle names of split bam file."""
    ## return the list of sub names
    return ['mapped','placed','unmapped',mito,'collisions'] 

## Ftn for commenting out command lines for debuging
def debuglines(intxt):
    ## Initilizse new lines 
    newlines = []
    ## Iterate thru the input txt lines 
    for l in intxt: 
        ## If the first chracter is already a comment like #SBATCH
        if l[0] == '#':
            newlines.append(l)
        ## If it is an echo statment, leave it as is 
        elif l.split(' ')[0]=='echo':
            newlines.append('sleep 60\n'+l)
        else: ## Othewise, comment out the lines 
            newlines.append('##'+l)
    ## Return the commented out lines 
    return newlines 

## Ftn to write to file
def writetofile(inpath,intxt,debug,mode='w'):
    """Opens a file to write lines to file."""
    ## Modify the input text lines if in debug mode
    intxt = debuglines(intxt) if debug else intxt
    with open(inpath,mode) as ofile:
        ofile.writelines(intxt)
    ## Return the path
    return inpath

## Write a ftn for formating lines for printing
def display(inlines,sep='\n'):
    """Formats input lines for printing to screen."""
    ## Joins the lines in limes
    return sep.join(inlines)

## Ftn for writing a set of read names
def writeset(outfile,inset):
    """Writes a file of a give set of read names."""
    ## Open a file for writing
    return writetofile(outfile, display(list(inset))+'\n', False)

## --------------------------------- SAMTOOLS VERSION ----------------------------------- ## 
## Ftn for calling the commands to shell
def submitcall(command):
    """Submits a command to the shell via subprocess."""
    ## Submit the input command
    return subprocess.call(command,shell=True)

## Write ftn for getting samtoosl version 
def getsamv(com='samtools --version',coding='utf-8'):
    """Gathers the version of samtools as a string."""
    ## Returns the output of samools version
    return subprocess.check_output(com,shell=True).decode(coding).split('\n')[0].split(' ')[-1]

## Write ftn for spliting and making an int
def splitint(x,sep='.'):
    """Splits an input string x on SEP and converts values to integers."""
    ## Return the array of split ints
    return [int(z) for z in x.split(sep)]

## Write the sam error 
sam_error = "ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again."

## Check that we have the correct version of samtools
def checksam(version='1.15.1'):
    """Returns a boolean, checking if the version of samtools is correct."""
    ## Return the bool of the minimum difference
    return min([ a-b for a,b in zip(splitint(getsamv()),splitint(version))])>=0

## ---------------------------------- PYSAM FUNCITONS ---------------------------------- ## 
## Ftn for checking if we have a bam file
def isbam(inbam):
    """Checks if the input file is a .bam file"""
    ## Return the bool
    return ((inbam.split('.')[-1] == 'bam') and isfile(inbam))

## Write a ftn for checking if index is there
def hasix(inbam):
    """Checks if the index file .bai of an input bam file exists."""
    ## Return the boolean
    return isfile(inbam+'.bai') or isfile(inbam+'.csi')

## Write ftns used in analysis ftn for loading a bam
def loadbam(inpath,bool=True):
    """Reads in bam file in binary format."""
    ## Return the alignment file
    return pysam.AlignmentFile(inpath,'rb',check_sq=bool)

## Ftn for gathering alignment recores from pybam object
def getrecs(pybam):
    """Gathers the records within an input python bam object."""
    ## returns the alignment recores from pybam object
    return [rec for rec in pybam] 

## Ftn for getting set of reads
def readset(pybam):
    """Returns the set of read names from an input python bam object."""
    ## return the set
    return set([rec.qname for rec in pybam]) 

## Gather the pairs and the records
def getpairs(records):
    """Returns the list of paired and singlton alignment records from an input list."""
    ## Initilize the pairs and singls
    pairs,singles = [], []
    ## Iterat thru the records
    for rec in records:
        if rec.is_paired: ## If it was paired append to pairs
            pairs.append(rec)
        else: ## Otherwise append the record to singlton
             singles.append(rec)
    ## Check our work, the lineths should match 
    assert lencomp(pairs + singles, records), 'ERROR: The number of paired and singlton records given as input do not match the output in splitting!\n'
    ## Return the pairs and singletons
    return pairs, singles

## Ftn for checking if a record is primary
def isprimary(pyrec):
    """Checks if a record is not a primary or secondary alignment."""
    ## Return if not supplementary or secondary
    return not (pyrec.is_supplementary or pyrec.is_secondary)

## Ftn for checking if a record is mapped
def ismapped(pyrec):
    """Returns a boolean on a py-alignment-record and its mate both being mapped."""
    ## return the combined boolean
    return (not pyrec.is_unmapped) and (not pyrec.mate_is_unmapped)

## Ftn for checking if a record is placed, mate is mapped, it is not mapped, 
def isplaced(pyrec):
    """Returns a boolean on a py-alignment-record if it is mapped but made is not."""
    ## Returns a boolean if a read is placed
    return ((not pyrec.is_unmapped) and pyrec.mate_is_unmapped) or (pyrec.is_unmapped and (not pyrec.mate_is_unmapped))

## Ftn for checking if an aligment record is a mtDNA
def ismtDNA(pyrec,mito):
    """Reaturns a boolean on a py-alignment-record or its mate being mapped to the mtDNA."""
    ## Return the boolean
    return (pyrec.reference_name == mito) or (pyrec.next_reference_name == mito)

## Ftn for checking if a record is unmapped
def unmapped(pyrec):
    """Returns a boolean on an unmapped py-alignment-record and its mate both being un-mapped."""
    ## Boolean if both the record and its mate are unmapped
    return pyrec.is_unmapped and pyrec.mate_is_unmapped

## Ftn for removing .txt and adding .bam
def txttobam(intxt):
    """Formats a bam file name from an input txt file name."""
    ## Return the formated bam name
    return intxt.split('.txt')[0] + '.bam'

## Ftn for counting 
def bamcount(inbam,countfile,threads):
    """Formats the view commadn to count primary aligments within a bam file."""
    ## Format the view command to count 
    return f'samtools view -c -F 256 -f 64 -@ {threads} {inbam} > {countfile}\n'

## Samtools view command to filter by read name 
def bambyreadname(inbam,readfile,threads,opts='-b'):
    """Submits a call to samtools view to filter input by read name."""
    ## format outbam
    outbam = txttobam(readfile)
    ## Formats the call to samtools for makign a bam file from list of read names and index said file
    return f'samtools view {opts} -N {readfile} -@ {threads} {inbam} | samtools sort -@ {threads} -o {outbam} -O BAM --write-index\n', outbam 

## Ftn for callign samblaster
def samblaster(inbam,outbam,report,threads):
    """Formats a call to samtools and samblaster given inputs."""
    ## Return the samtools and samblaster command 
    return f'samtools sort -@ {threads} -n {inbam} | samtools view -@ {threads} -Sh - -O SAM | samblaster --ignoreUnmated -M 2>> {report} | samtools view -@ {threads} -Shb | samtools sort -@ {threads} - -o {outbam} -O BAM --write-index\n'

## Ftn for formating flag
def dupflag(bool):
    """returns a lower or upper case F for samtools view commands given boolean."""
    ## Format the samtool flag
    return '-f' if bool else '-F'

## Ftn for formating sort command removing duplicates
def dupsort(bamin,bamout,threads,keep=False,dupint=1024):
    """Formats samtoosl view and sort command to keep or remove duplicates given flags in opts."""
    ## Return the formated txt 
    return f'samtools view -@ {threads} {dupflag(keep)} {dupint} -Shb {bamin} | samtools sort -@ {threads} - -o {bamout} -O BAM --write-index\n'

## Ftn for filtering primary
def getprimary(inbam,mapq,threads,outbam):
    """Formats a samtools view command with the -e flag to return primary aligments from an input bam file (inbam) at a given mapping quality (mapq) and returns aligments within an output bam file (outbam)."""
    ## Return the samtools commands 
    return f'samtools view -@ {threads} -f 3 -F 1024 -b {inbam} -e "!([XA] || [SA]) && mapq>={mapq}" | samtools sort -@ {threads} - -o {outbam} -O BAM --write-index\n'

## ------------------------------------------------------- HIC FUNCTIONS ---------------------------------------------------------------------- ## 
## Ftn for fetching cigar string
def getcigar(pyrec):
    """Returns a cigar tuple based on read orientation in pyrec."""
    ## return the cigar str based on strand orientation
    return pyrec.cigartuples[::-1] if pyrec.is_reverse else pyrec.cigartuples[:]

## Reurn the match sumation of py records
def getmatchsum(cigartuples):
    """Returns the summed match score across cigar tuples from alignment record."""
    ## Set match sum to start at zero
    match_sum = 0
    ## Iterate over the cigar tuple
    for m, c in cigartuples[:]:
        ## If m is a match to the cmatch from pysam, break from loop
        if (m == pysam.CMATCH):
            break
        ## Add to match sum 
        match_sum += c
    ## Return the match sum
    return match_sum

## Ftn for parsing the correct split of an input record
def correctsplit(pyrecs):
    """Calculates a match position score given a split alignment record."""
    mapped_correctly = []
    for rec in pyrecs:
        cigar = getcigar(rec)
        match_summed = getmatchsum(cigar)
        mapped_correctly.append(match_summed)
    min_rec_idx = mapped_correctly.index(min(mapped_correctly))
    ## Return the record with the minimum matching positional index
    return pyrecs[min_rec_idx]

## Ftn for getting correct alignment
def getalignment(pyrecs):
    """Returns a primary alignemnt or the correct split alignment."""
    ## Return a correct split or the primary alignments
    return correctsplit(pyrecs) if (len(pyrecs)>1) else pyrecs[0]

## Ftn for formating strand from pyrecord
def strand(orientation):
    """Returns a -1 or 1 if the orientation is reversed (True) or not (False)."""
    ## return the orientation 
    return -1 if orientation else 1

## Ftn for orienting a pair of reads
def orderpair(recA,recB):
    """Orders input records A and B by reference position."""
    ## If the reads on the same chormosome 
    if (recA.reference_id == recB.reference_id):
        ## If the position of A is less than or equal to B
        if recA.pos <= recB.pos:
            rec1, rec2 = recA, recB
        else: ## otherwise make B the first
            rec1, rec2 = recB, recA 
    else: ## If the reads are on different chromosomes 
        if recA.reference_id < recB.reference_id: 
            rec1, rec2 = recA, recB
        else: ## Otherwise make record B the first 
            rec1, rec2 = recB, recA 
    ## Return the ordered recs
    return rec1, rec2 

## Ftn for formating pyrec(s) in a line seen in pair-end bed file
def bamtobedpe(pyrecA,pyrecB):
    """Gathers and formats mapping information for an input pair of reads stored as py-records from py-bam object."""
    ## Order the records, left to right. 
    pyrec1, pyrec2 = orderpair(pyrecA,pyrecB)
    ## Gather the alignment information for input pair of reads
    qname            = pyrec1.qname                                           ## set the read name
    mapq1, mapq2     = pyrec1.mapq, pyrec2.mapq                               ## Mapping quality
    chr1, chr2       = pyrec1.reference_name, pyrec2.reference_name           ## Reference name
    chrn1, chrn2     = pyrec1.reference_id, pyrec2.reference_id               ## Reference id 
    left1, left2     = pyrec1.pos, pyrec2.pos                                 ## chromosome position, start from left
    right1, right2   = left1 + pyrec1.qlen, left2+pyrec2.qlen                 ## end position (adding read length to start)
    strand1, strand2 = strand(pyrec1.is_reverse), strand(pyrec2.is_reverse)   ## Format strand 
    ## Format and return the row
    return f'{chr1} {chrn1} {left1} {right1} {mapq1} {strand1} {chr2} {chrn2} {left2} {right2} {mapq2} {strand2} {qname}'
## End of file 