#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## ----------------------------- v 8.0.0 ------------------------------ ## 
##         PYBWATOOLS: Python functions for working with BWA 
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
## ----------------------------------- FORMATING FUNCTIONS --------------------------------- ##
## Ftn for echo bwa finish to a log
def bwaecho(o,l=None):
    """Formats an echo command to print finishing statment of bwa alignment to log."""
    ## Format the message given the output (o) file name and log (l)
    return f'echo Finished alignment and indexing of split: {o} > {l}' if l else ''

## Ftn for mapping reads in hic-mode as suggested by phase genomics and juicer
def bwamem_hic(r1,r2,ref,outbam,log,vmode=1,threads=4,opts='-5SPM'):
    """Formats a bwa mem command for paired-end hi-c sequencing reads."""
    ## Set the bwa mem option in hic mode and format the bwa mem call, submit to shell via submit command
    return [f'bwa mem {opts} -v {vmode} -t {threads} {ref} {r1} {r2} 2>> {log} | samtools view -Shb -@ {threads} | samtools sort -@ {threads} - -o {outbam} -O BAM --write-index\n', bwaecho(outbam,log)]

## Ftn for formating a bwa mem for a single read file
def bwamem_single(read,ref,outbam,log,vmode=1,threads=4):
    """Formats a bwa mem command for a single read alignment given the read and reference name."""
    ## Format the bwa meme command for singleton reads
    return [f'bwa mem -v {vmode} -t {threads} {ref} {read} 2>> {log} | samtools view -Shb -@ {threads} | samtools sort -@ {threads} - -o {outbam} -O BAM --write-index\n', bwaecho(outbam,log)]

## Write ftn for formating bwa mem command with paired reads
def bwamem_paired(r1,r2,ref,outbam,log,vmode=1,threads=4,opts='-M'):
    """Formats a bwa command and requires as input paired-reads, a reference name, an output bam file, a log file, experiment type."""
    ## Set the bwa mem option based on experiment type and format the bwa mem call, submit to shell via submit command
    return [f'bwa mem {opts} -v {vmode} -t {threads} {ref} {r1} {r2} 2>> {log} | samtools view -Shb -@ {threads} | samtools sort -@ {threads} - -o {outbam} -O BAM --write-index\n', bwaecho(outbam,log)]
## End of file 