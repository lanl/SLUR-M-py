#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
#######################################
##     SLURPY PIPELINE (v 8.0.0)     ##
#######################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## -------------------------------------------------------------- MODULE LOADING ----------------------------------------------------------------- ## 
## Load in glob
from glob import glob 

## Bring in sub-process mod
import subprocess

## Bring in basename, get file size,  and pathexists 
from os.path import basename, exists as pathexists, getsize as getfilesize

## Bring in make dirs
from os import makedirs, remove 

## Bring in our costum ftns from pysamtools 
from pysamtools import ifprint, splitbam, writetofile, makelist, checksam, dictzip, bambyreadname, display, samblaster, getprimary, txttobam, outnames

## Load in bwa ftns
from pybwatools import bwamem_hic, bwamem_paired, bwamem_single

## ------------------------------------------------------- DEFAULT VARIABLE SETTING ---------------------------------------------------------------- ## 
## Set the version and description of slurpy pipeline 
myversion = '8.0.0'
slurpy_descr = f'Calls the slurpy pipeline ( v {myversion} ) for alignment of pair-end reads to a reference genome.\nExample call:\n\t./slurpy -e atac -r ../../REFERENCES/ENCODEREF/AUTOSOMES/GRCh38.fasta'

## Set directory names 
debugdir   = 'debug'        ## Hold logs for debug 
fastqdir   = 'fastqs'       ## The directory holding fastq files 
aligndir   = 'aligned'      ## Holds final aligments 
splitsdir  = 'splits'       ## Temporary dir for split fastq 
comsdir    = 'commands'     ## Folder for holding all command files 
macs2dir   = 'macs2'        ## Has results from macs2 
hicdir     = 'hic'          ## Has hic resluts 
diagdir    = 'diagnostics'  ## Plots for diagnostics are held here 
bamtmpdir  = 'bamtmp'       ## A temporary dir for hodling bam files from split fastq aligments 
scriptsdir = './SLURPY'     ## The script directory holding this file

## Group the dirs 
grouped_dirs = [debugdir,aligndir,splitsdir,comsdir,diagdir,bamtmpdir]

## ----------------------------------------------------- SLURPY FUNCTIONS  ---------------------------------------------------------------- ## 
## Ftn for returning a sorted glob
def sortglob(wildcard):
    """Retuns the sorted glob of input wild card."""
    ## Returns the sorted glob of input
    return sorted(glob(wildcard))

## Ftn for getting sample name
def getsamplename(read1,spliton='_R1_'):
    """Get sample name from a fastq.gz file."""
    ## Return the name 
    return basename(read1).split(spliton)[0]

## Ftn for returning just the bam file name without the .bam
def basenobam(inbam):
    """Returns the basename of an input bam file without the .bam file extension."""
    ## Return the split basename
    return basename(splitbam(inbam))

## Ftn for takeing basne anem 
def basenosh(insh):
    ## Returns the basename of in .sh file with no .sh
    return basename(insh).split('.sh')[0]

## Write a ftn to check the file size meets our threshold
def checkfilesize(filepath,threshold=5):
    """Returns a boolean on the input file path checking if a given file is larger (in bytes) than threshold."""
    ## Return the boolean of file size
    return getfilesize(filepath) > threshold

## Write ftn for checking file is a real deal!
def fileexists(filepath):
    """Checks the existance and size of an input file (foodstuffs)."""
    ## Via os, check if the path to "food" exists and meets our size threshold
    return pathexists(filepath) and checkfilesize(filepath)

def addfslash(inpath):
    """Corrects input inpath with a forward slash."""
    ## Returns a path if the last chracter is a forward path 
    return inpath if (inpath[-1] == '/') else inpath + '/'

## Write ftn for retreving the first read in pair
def getread1(fpath,sep='_R1_',fend='.fastq.gz'):
    """Lists the first reads (in fastq.gz format) in pair given the fastq directory."""
    ## Return the frist fastq in pair
    return sortglob(f'{fpath}*{sep}*{fend}')

## Ftn for checking pairs of fastq files
def ispared(r1,r2):
    """Checks that both r1 and r2 are paired strings and not None type."""
    ## Return the paired fastq files if they have the same type and are not None type.
    return [(r1,r2) for (r1,r2) in zip(r1,r2) if (type(r1)==type(r2)) and type(r2)]

## Ftn for formating read 2
def setread2(r,seps):
    """Format read 2 from read 1 in r on seporators in seps."""
    ## Return the formated read 2
    return seps[1].join(r.split(seps[0]))

## Write ftn for pairing fastqs
def getfastqs(fpath,seps=('_R1_','_R2_'),fend='.fastq.gz'):
    """Lists fastq files in pair from input path."""
    ## Gather read 1 fastq
    read1 = getread1(fpath, sep=seps[0], fend=fend)
    ## Check the pairs and get a paired list of reads
    pairs = ispared(read1, [setread2(r,seps) if fileexists(setread2(r,seps)) else None for r in read1])
    ## Return the pairs and the message of warn
    return pairs

## Ftn for calling sub-process
def submitsbatch(pathtoscript):
    """Checks output from sub-process submission of input path-to-script assuming it is an sbatch command. Returns the job id"""
    ## call sub process
    return subprocess.check_output(pathtoscript, shell=True).decode("utf-8").split('\n')[0].split(' ')[-1]

## Ftn for making items within a list into strings
def makestrs(inlist):
    """Formats all items within input list into a string."""
    ## Makes all items within input list into a string
    return [str(i) for i in inlist]

## Formats digit 
def formatn(n):
    """Formats input digit to have 4 decimal places."""
    ## Format an integer to 4 deciimal places
    return "{:04d}".format(n)

## Ftn for checking if a file is bwa indexed
def isbwaix(inref,indexends = ['amb', 'ann', 'bwt', 'pac', 'sa']):
    """Checks if the the input reference file is indexed."""
    ## Check to make sure the reference exists
    assert fileexists(inref), "ERROR: Unable to locate input reference file on path: %s"%inref
    ## Check if each of these files exist: amb, ann, bwt, pac and sa
    return sum([pathexists(inref+'.'+fe) for fe in indexends]) == len(indexends)

## Ftn for formating an sbatch text
def sbatch(nameojob,cpus,cwd,partition=None,nodes=1,tasks=1,runtime='200:00:00'):
    ## return the formated sbatch txt
    settings = ['#!/usr/bin/env bash\n',
                '#SBATCH --job-name=%s\n'%basenosh(nameojob) if nameojob else '\n',
               f'#SBATCH --output={debugdir}/%x.%j.out' +'\n',
               f'#SBATCH --error={debugdir}/%x.%j.err' + '\n',
                '#SBATCH --nodes=%s\n'%str(nodes),
                '#SBATCH --ntasks-per-node=%s\n'%str(tasks),
                '#SBATCH --cpus-per-task=%s\n'%str(cpus),
                '#SBATCH --time=%s\n'%runtime]
    ## Add the current working dir 
    settings = settings + ['#SBATCH --chdir=%s\n'%cwd] if cwd else settings
    ## Add the partition 
    settings = settings + ['#SBATCH --partition=%s\n'%partition] if partition else settings
    ## return settings 
    return settings

## Ftn for the begining of a call to fastp
def fastcut(i,I,o,O,r,n):
    """Formats the start of a call to fastp given inputs (i,I), outputs (o,O), and the report name (r)."""
    ## Format the fall and return
    return f'fastp -i {i} -I {I} -o {o} -O {O} -j {r}.{n}.json -h {r}.{n}.html'

## Ftn for returning a boolean if in hic
def inhic(inexp):
    """Boolean if the input experiment is hi-c."""
    ## Return the hic 
    return inexp == 'hic'

## Write ftn for logging fastp has finished
def fastdry(r1,r2,report):
    """Formats an echo command for logging the completion of a fastp split."""
    ## reformat the report 
    report = report if report.split('.')[-1] == 'txt' else report + '.txt'
    ## Format and return the echo command
    return f'echo Finished filtering and splitting on: {r1} {r2} >> {report}' + '\n'

## Ftn for formating fastp command to filter and split reads
def fastpeel(r1,r2,exp,w,s,z=4,options=['--dont_eval_duplication','--disable_length_filtering','--disable_adapter_trimming','--disable_quality_filtering','--disable_trim_poly_g'],script='fastp',toremoveend='.toberemoved.fastq.gz',singleend='.singletons.fastq.gz',failend='.failed.fastq.gz'):
    """Formats calls to fastp given input fastq files."""
    ## Gather r1 and r2 basename
    r1_bn, r2_bn     = basename(r1).split('.fastq')[0],    basename(r2).split('.fastq')[0]
    ## Format the temp splits to be remove
    temp1, temp2     = f'{splitsdir}/{r1_bn}{toremoveend}', f'{splitsdir}/{r2_bn}{toremoveend}'
    ## Format the single read files
    single1, single2 = f'{splitsdir}/{r1_bn}{singleend}',   f'{splitsdir}/{r2_bn}{singleend}'
    ## Format the splits
    split1, split2   = f'{splitsdir}/{basename(r1)}',       f'{splitsdir}/{basename(r2)}'  
    ## Parse the report, failed name, and experiment 
    report, failed = reportname(getsamplename(r1),script),  f'{splitsdir}/{getsamplename(r1)}{failend}'
    ## Format the initial run of fastp to filter reads and trim addaptors
    wash = fastcut(r1,r2,temp1,temp2,report,0) + f' --unpaired1 {single1} --unpaired2 {single2} --failed_out {failed} --thread {w} ' + ' '.join(options[:2]) + '\n'
    ## Call fastp again to split the files
    cut = fastcut(temp1,temp2,split1,split2,report,1) + f' -s {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'
    ## Call fastp again to split the files
    slice = fastcut(r1,r2,split1,split2,report,1) + f' -s {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'
    ## Call the fast dry command and format the remove command 
    dry, throwout = fastdry(r1,r2,report), f'rm {temp1} {temp2}\n'
    ## Format and return the command based on the experiment type (ie if it is hi-c or not)
    return [slice,dry] if inhic(exp) else [wash,cut,throwout,dry], report 

## Ftn for formating report
def reportname(inbam,script):
    """Formats the report name given the input bam file name and type of script (str)."""
    ## Return the report name
    return f'{debugdir}/{basenobam(inbam)}.{script}.log.txt'

## Ftn for formating single epfastq files for bwa mem in hic mode
def prepbwamem(r1,r2,refix,t,mode):
    """Formats a bwa mem command given the paired reads, reference index, thread count, and mode."""
    ## Format the bam file name and report name
    bam, report = f'{bamtmpdir}/{getsamplename(r1)}.bam', f'{debugdir}/{getsamplename(r1)}.bwa.log.txt'
    ## Format the bwa mem command 
    if inhic(mode): ## Hi-c analysis mode 
        makebams = bwamem_hic(r1,r2,refix,bam,report,threads=t)
    ## If thie mdoe is single tons 
    elif (mode == 'singletons') and (r2 is None): 
        makebams = bwamem_single(r1,r2,refix,bam,report,threads=t)
    ## We are in pair end mode 
    else: ## We are in paired end mode 
        makebams = bwamem_paired(r1,r2,refix,bam,report,threads=t)
    ## Return the list of bams to make and the count of those made
    return bam, makebams, report  

## Ftn for formating samblaster command
def markduplicates(inbam,threads,script='mark'):
    """Formats and submits a samtools and samblaster command to mark duplicates on input bam file (i) and saves to output (o)."""
    ## Format the output bams and the report name 
    outbam, report = f'{aligndir}/{basenobam(inbam)}.marked.bam', reportname(inbam,script) 
    ## Format the sam-blaster and echo command
    blast_command, echo_command = samblaster(inbam,outbam,report,threads), f'echo Finished marking duplicates in {outbam} >> {report}\n'
    ## Return the samblaster command and ecco chommand 
    return outbam, [blast_command, echo_command], report

## Ftn for formating echo split
def splitecho(inbam,report,script):
    """Formats an echo statment for the split command."""
    ## Return the formated command 
    return f'echo Finished splitting {inbam} using {script} >> {report}\n'

## Ftn for formating call to our split command
def splitcommand(bam,mito,exp,threads,splitsize=10000,script='split'):
    """Formats and returns a command to our protocol for splitting an input bam file."""
    ## Call the report and echo command
    report = reportname(bam,script)
    ## Set command based on experiment type
    if inhic(exp): ## Format the hic splitting command
        commands = [f'{scriptsdir}/hicspliter.py -b {bam} -M {mito} >> {report}\n', splitecho(bam,report,'hicspliter.py')]
    else: ## Otherwise format the bam splitting script 
        ## Set the output names fromt he bam file 
        mapped_txt, placed_txt, mito_txt, unmap_txt, bedpe_txt = outnames(bam,mito)
        ## Format the command for splitting the input bam file, index the output bam file, and samtools commands to form bam files from read names
        commands = [f'{scriptsdir}/readspliter.py -b {bam} -M {mito} >> {report}\n'] + [bambyreadname(bam,r,threads)[0] for r in [mapped_txt, placed_txt, mito_txt]] + [splitecho(bam,report,'readspliter.py')]
    ## Return the command
    return commands, report 

## Write ftn for making a directory
def dirmaker(dirpath):
    """Makes a directory given an input path."""
    ## Return the os command
    return makedirs(dirpath,exist_ok=True)

## Write ftn for making folders used here
def makedirectories(folderpaths):
    """Generates the directories listed in folderpaths."""
    ## Iteratively make the dir in the folder paths
    return [dirmaker(f) for f in makelist(folderpaths)]

## Write ftn for making head path
def headpath(inpath):
    """Returns the head path of a path of type string."""
    ## Split the path take up to the last item and reformat
    return '/'.join(inpath.split('/')[:-1]) 

## Ftn for formating bwa index
def bwaindex(refpath,script='index'):
    """Formats command to index a reference via bwa."""
    ## return the index commands
    return [f'bwa index {refpath}\n', f'echo Finished indexing reference on path {refpath} > {reportname(refpath,script)}\n'], reportname(refpath,script)

## Ftn for formating job ids
def formatids(cdf,op,joinon=','):
    """Formats job ids within the command dataframe (cdf) given an operation."""
    ## Return the joined list of ides 
    return joinon.join([str(k) for k in cdf[(cdf.Operation.isin(op))].JobID.tolist() if len(str(k)) > 0])

## Ftn for fromating dependencies
def formatdepends(jobids,depend='afterok'):
    """Formats job ids and dependencys for an sbatch submission."""
    ## Return the formated text
    return f' --dependency={depend}:{jobids} ' if (len(jobids) > 0) else ' '

## Ftn for formating job name 
def jobname(args):
    """Formats an input job name, joining inputs in args with a period."""
    ## Return the job name
    return '--job-name=%s'%'.'.join([str(a) for a in args])

## Ftn for submission of jobs with dependencies 
def submitdependency(command_df,operation,dependent,timestamp,clusterpart,bylast=False,group='Sample'):
    """Formats and submits sbatch jobs from command dataframe based on operations and dependents"""
    ## Initilzse cap and lists
    subsbatchs, dependent = [], makelist(dependent)
    ## group the 
    for sample_ix, sample_df in command_df[(command_df.Torun==0)].groupby(group):
        ## Gather the operaiton and dependent ix
        operations, dependents = sample_df[(sample_df.Operation==operation)], sample_df[(sample_df.Operation.isin(dependent))]
        ## Gather the operations indexes 
        op_ixs = operations.index.values
        ## Iterate thru the op indexs
        for m in op_ixs:
            ## gather the dependent ix
            dependent_ix = [i for i in dependents[(dependents.index < m)].index.tolist()]
            ## Take the last index if we are doing it by lst
            dependent_ix = dependent_ix[-1:] if bylast else dependent_ix
            ## Set the row
            row = command_df.loc[m,:]
            ## Format the bwa ids
            djobids, jobparts = formatids(command_df[(command_df.index.isin(dependent_ix))],dependent), jobname([operation,m,clusterpart,timestamp])
            ## Set the sbatch text
            intext = f'sbatch{formatdepends(djobids)}--partition={clusterpart} {jobparts} {row.Jobfile}'
            ## print the comamnd
            ifprint(intext,debug)
            ## Append to sub list
            subsbatchs.append(intext)
            ## Submit the job and append the job id to the dataframe 
            command_df.loc[m,'JobID'] = fakejobid + m if (debug and runlocal) else submitsbatch(intext)
    ## Return submitted sbatchs 
    return subsbatchs

## Ftn for formating merging of bam files
def mergebam(bam,wildcard,threads,script='merge'):
    """Formats a samtools merge command."""
    ## Format path of out put bam file
    outbam = f'{aligndir}/{bam}'
    ## Format report name and the merge-bam command
    report, merge_bam_command = reportname(outbam,script), f'samtools merge -f -@ {threads} -o {outbam} --write-index {bamtmpdir}/{wildcard}\n'
    ## Format the echo command and count command 
    echo_merge_command = f'echo Finished merging bam files into {outbam} >> {report}\n'
    ## Return the formated merge command
    return [merge_bam_command,echo_merge_command], report

## Ftn for concatonating
def concatcom(outfile,inwildcard,script='merge'):
    """Formats a command to merge file."""
    ## format the report
    report = reportname(outfile,script)
    ## Return the formated commands
    return [f'cat {inwildcard} > {outfile}\n',f'echo Finished concatenating files into {outfile} >> {report}\n'], report 

## Ftn for filtering bam fie 
def filterbam(inbam,mapq,threads,script='filter'):
    """Formats a command to filter an input bam file seperating aligments on quality."""
    ## Format the report name and out bam name 
    report, outbam = reportname(inbam,script), splitbam(inbam) + f'.primary.q{mapq}.bam'
    ## Format filter command and the echo command 
    bam_filter_command, echo_command = getprimary(inbam,mapq,threads,outbam), f'echo Finished filtering {inbam} at mapq of {mapq} >> {report}\n'
    ## Format and return commands
    return outbam, [bam_filter_command, echo_command], report

## Ftn for setting genomesize
def genomesize(inputsize,referencepath,mtDNA,sep='\t',header=None):
    """Calculate the genome size given an input size or path to reference genome .fai file."""
    ## set the genome size, load in the ref fai file and gahter the chromsome lengths 
    if inputsize: ## If the input size is none 
        genome_size = inputsize
    else: ## Load in the size dataframe 
        size_df = pd.read_csv(referencepath+'.fai',sep=sep,header=header)
        ## Calculate genome size
        genome_size = size_df[~(size_df[0].isin(makelist(mtDNA)))][1].sum()
    ## Return the size
    return genome_size

## Ftn for formating command control dataframe and restarting
def commandcontrol(commands,toreset,pipelinesteps,rerunfrom,bwaid,cols=['Jobfile','Dependency','Sample','Experiment','Operation','AfterID','Report']):
    """Generates a command and control dataframe for running pipeline."""
    ## Make into a dataframe
    commanddf = pd.DataFrame(commands,columns = cols)
    ## We dont need the dependency column or after id (remove these in next version )
    commanddf.drop(['Dependency','AfterID'],axis=1,inplace=True)
    ## Add to run column and job id 
    commanddf['Torun'], commanddf['JobID'] = 0, ''
    ## If a hard reset was called, remove the previous file reports 
    if_hard_reset = [(remove(fr) if fileexists(fr) else None) for fr in commanddf.Report.tolist()] if toreset else None 
    ## Check if the reports exist, iterate of the rows of command df
    for rix,row in commanddf.iterrows():
        commanddf.loc[rix,'Torun'] = 1 if fileexists(row.Report) else 0
    ## For completeness, set the bwa index to run to one 
    commanddf.loc[(commanddf.Operation=='bwaix'),'Torun'], commanddf.loc[(commanddf.Operation=='bwaix'),'JobID'] = 1, bwaid 
    ## If we are re running any part of the pipeline we will re code from here
    if rerunfrom: ## Set the dict 
        rerun_dict = dictzip(pipelinesteps,[pipelinesteps[i:] for i in range(len(pipelinesteps))])
        ## Recode the command dataframe to run 
        for r in rerun_dict[rerunfrom]:
            commanddf.loc[(commanddf.Operation==r),'Torun'] = 0
    ## Return the command dataframe and hard reset 
    return commanddf, if_hard_reset

## ------------------------------------------------------- DEFAULT VARIABLE SETTING ---------------------------------------------------------------- ## 
## Set default split size
splitsize   = 16
bwathreads  = 4
samthreads  = 4
parallelbwa = splitsize
part        = 'tb'
map_q_thres = 30

## Set the mitocondrial contig
mito_contig = 'chrM'

## Set protocols and pipe lien steps
protocols        = ['atac','chip','wgs','hic']
joined_protocols = display(protocols,sep=', ')
pipeline_steps   = ['fastp', 'bwa', 'split', 'merge', 'mark', 'filter','macs2','count','clean']

## ----------------------------------------------------------- HELP MESSAGES ------------------------------------------------------------------------ ## 
## Set help messages
e_help = "The experiment type (default: wgs). Valid options include: %s."%joined_protocols
r_help = "Path to input reference bwa index used in analysis." 
F_help = "The number of splits to make for each pair of input fastq files (default: %s). Controls the total number of splits across the run."%splitsize
f_help = "The number of threads used in fastp to split input fastq files (default: %s). Note: must be an even multiple of the number of splits."%splitsize
b_help = "The number of threads used per bwa alignment on split input fastq files (default: %s)."%bwathreads
n_help = "Run name used to name output files. Default behavior is to use the current parent directory."
M_help = "Name of the mitochondrial contig (default: %s)."%mito_contig
B_help = "Number of parallel bwa alignments to run (default: %s). Controls the number of bwa jobs submitted at once to slurm."%parallelbwa
P_help = "The type of partition jobs formatted by slurpy run on (default: %s)."%part
R_help = "Step within the pipeline to re-run from. Options (in order) include: %s."%', '.join(pipeline_steps)
Q_help = "Mapping quality threshold to filter alignments (default: %s)."%map_q_thres
c_help = "Path to control or input bam files used in ChIP-seq experiments."
g_help = "Size of the genome being analyzed, used as parameter for macs2. Inputs can be integers in bp or two letter short hand, for e.g. hs for homo sapiens."

## Set help messages for bollean vars
restart_help  = "Flag to force the pipeline to reset and run from start."
runlocal_help = "Disables sbatch submission and submits the script via bash to a local os."
debug_help    = "A flag to run in verbose mode, printing sbatch commands. Default behavior is false."
mark_help     = "Pass this flag to skip marking and removing duplicates. Default behavior is false (conduct duplicate marking)."
broad_help    = "Flag to call broad peaks using the --broad-cutoff=0.1 setting in macs2. See macs2 callpeak --help for more details."
clean_help    = "If included will run clean up script at end of run. The default behavior is false, can be run after pipeline."

## Set metavars
c_metavar = './path/to/control.bam'
g_metavar = 'bp'

## Write error message
not_sam_err = "ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again."

## ------------------------------------------------- MAIN SCRIPT & ARGUMENT PARSING --------------------------------------------------------------- ## 
## If the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = slurpy_descr)

    ## Add the required argument
    parser.add_argument("-r", "--refix",         dest="r", type=str,  required=True, help = r_help,  metavar = './path/to/reference.bwaix'        )
    
    ## Add the default arguments
    parser.add_argument("-F", "--fastp-splits",  dest="F", type=int,  required=False, help = F_help, metavar = splitsize,   default = splitsize   )
    parser.add_argument("-B", "--parallel-bwa",  dest="B", type=int,  required=False, help = B_help, metavar = parallelbwa, default = parallelbwa )
    parser.add_argument("-P", "--partition",     dest="P", type=str,  required=False, help = P_help, metavar = part,        default = part        ) 
    parser.add_argument("-M", "--mtDNA",         dest="M", type=str,  required=False, help = M_help, metavar = mito_contig, default = mito_contig )
    parser.add_argument("-Q", "--map-threshold", dest="Q", type=int,  required=False, help = Q_help, metavar = map_q_thres, default = map_q_thres )

    ## Set number of threads 
    parser.add_argument("-f", "--fastp-threads", dest="f", type=int,  required=False, help = f_help, metavar = splitsize,   default = splitsize   )
    parser.add_argument("-b", "--bwa-threads",   dest="b", type=int,  required=False, help = b_help, metavar = bwathreads,  default = bwathreads  )
    parser.add_argument("-s", "--sam-threads",   dest="s", type=int,  required=False, help = b_help, metavar = samthreads,  default = samthreads  )

    ## Set ftns for macs2
    parser.add_argument("-n", "--run-name",      dest="n", type=str,  required=False, help = n_help, metavar = 'name',      default = None        )
    parser.add_argument("-c", "--controls",      dest="c", nargs='+', required=False, help = c_help, metavar = c_metavar,   default = None        )
    parser.add_argument("-g", "--genome-size",   dest="g", type=str,  required=False, help = g_help, metavar = g_metavar,   default = None        )
    parser.add_argument("-e", "--experiment",    dest="e", type=str,  required=False, help = e_help, metavar = 'type',      default = 'wgs'       )

    ## Restart argument
    parser.add_argument("-R", "--rerun-from",    dest="R", type=str,  required=False, help = R_help, metavar = 'step',      default = None        )

    ## Set boolean flags 
    parser.add_argument("--restart",             dest="start",  help = restart_help,  action='store_true')
    parser.add_argument("--runlocal",            dest="runloc", help = runlocal_help, action='store_true')
    parser.add_argument("--debug",               dest="debug",  help = debug_help,    action='store_true')
    parser.add_argument("--skip-dedup",          dest="mark",   help = mark_help,     action='store_true')
    parser.add_argument("--broad",               dest="broad",  help = broad_help,    action='store_true')
    parser.add_argument("--clean",               dest="clean",  help = clean_help,    action='store_true')

    ## Set the paresed values as inputs
    inputs = parser.parse_args()
    ## ------------------------------------------------- PASS ARGUMENTS and SET VARIABLES ----------------------------------------------------------------- ## 
    ## Set variables
    reference_path = inputs.r       ## Set path to the reference genome 
    fastp_splits   = inputs.F       ## Number of splits in fastp 
    fastp_threads  = inputs.f       ## Number of fastp threads
    bwa_threads    = inputs.b       ## Number of threads in bwa alignments
    sam_threads    = inputs.s       ## Number of samtools threads 
    mito           = inputs.M       ## Set the mito contig name 
    bwa_runs       = inputs.B       ## Set the number of parallel runs of bwa 
    partition      = inputs.P       ## Set the partition 
    rerun          = inputs.R       ## Setp to rerun pipeline from 
    mapq           = inputs.Q       ## Set the mapping quality threshold 
    run_name       = inputs.n       ## The name of the samples 
    chip_control   = inputs.c       ## Set the input control for chip experimetn
    g_size         = inputs.g       ## Size of the genome   
    experi_mode    = inputs.e       ## The type of sequencing experiment 
        
    ## Set boolean vars 
    hardreset      = inputs.start   ## Resetart the slurpy run, removing previous
    runlocal       = inputs.runloc  ## Run script locaclly
    debug          = inputs.debug   ## Run in debug mode 
    skipduplicates = inputs.mark    ## Boolean to mark duplicates 
    ifbroad        = inputs.broad   ## Boolean to activate broader peak calling in macs2 
    ifclean        = inputs.clean   ## Flag to run clean up script 

    ## Hardset variables
    fakejobid = 400000             ## Fake job id used for debugging only 
    t2t_refpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCA_009914755.4_T2T-CHM13v2.0_genomic.named.fasta' ## Path to human t2t ref on canopus 

    ## ---------------------------------------- Cullen Roth Specific Substitutions ------------------------------------------------------- ## 
    ## Reset reference with presets if we are running as Cullen Roth (These can be edited if you are not me :-))
    reference_path = t2t_refpath if reference_path in ['T2T','t2t','T2t','t2T'] else reference_path

    ## ----------------------------------------------------- INITILIZATION ----------------------------------------------------------------------- ## 
    ## Check that the fastq path exists
    assert pathexists(fastqdir), "ERROR: Unable to detect a fastqs directory!"
    
    ## Check the versions of samtools, the user email is an email and the experiment mode is one we know
    assert checksam(), not_sam_err 
    #assert checkuser(user_email), f'ERROR: The input user email -- {user_email} -- is not a properly formated email address!'
    assert experi_mode in protocols, f'ERROR: The set experiment mode -- {experi_mode} -- is not in our list of portocols: {joined_protocols}'

    ## If needed reset that the fastp threads and splits such that they are a multiple of each
    if (fastp_splits%fastp_threads != 0):
        print(f'WARNING: The number of splits set for fastp ({fastp_splits}) is not an even multiple of the thread count ({fastp_threads})')
        ## Reset the thread count to half the split number
        fastp_threads = int(fastp_splits/2)
        ## Print we are resetting
        print(f'WARNING: The number fastp threads to conduct {fastp_splits} splits was reset to {fastp_threads}.')

    ## Check that if chip was called there is an input control
    if (experi_mode == 'chip'):
        assert (type(chip_control) == list) and (type(chip_control[0]) == str), f'ERROR: Inputs for chip control are not a type we recognize (should be a list of strings)!'
        ## Iterate over the controls inputs 
        for chip_con in chip_control:
            ## Check they exist 
            assert fileexists(chip_con), f'ERROR: Bad path to bam file! Unable to locate input control for chip experiment: {chip_con}'
    else: ## Otherwise do nothing 
        pass 

    ## Reformat clean boolean if clean was passed from restart
    ifclean = True if (rerun == 'clean') else ifclean

    ## --------------------------------------------------------- MODULE LOADING ------------------------------------------------------------------------------ ## 
    ## Load in the list dir ftn 
    from os import getcwd as gcwd

    ## Load pandas 
    import pandas as pd, time 

    ## Set the unique time stamp
    stamp = round(time.time())

    ## Load in macs2 ftns
    from pymacs2 import peakattack

    ## ------------------------------------------- DIRECTORY MAKING & TIME STAMP SUBMISSION ------------------------------------------------------------------- ## 
    ## Get the current working dir
    the_cwd = gcwd()

    ## If the run name is none
    run_name = run_name if run_name else basename(the_cwd)

    ## Set the broad pkeack
    broadpeak = '--broad' if ifbroad else ''

    ## Make the directories 
    makedirectories(grouped_dirs)

    ## If we are in atac or chip experiment mode
    if experi_mode in ['atac','chip']:
        makedirectories([macs2dir])
    else:
        pass 

    ## Gather and remove the old command files, error, and output logs
    [remove(f) for f in sortglob(f'{comsdir}/*.sh') + sortglob(f'{debugdir}/*.err') + sortglob(f'{debugdir}/*.out')]
    
    ## Generate the bwa indx command
    bwa_ix_jobname = f'{comsdir}/bwaindex.sh'

    ## Call the bwa command and its report 
    bwaix_coms, bwaix_report = bwaindex(reference_path)

    ## Write the bwa index command and sbatch to file 
    writetofile(bwa_ix_jobname, sbatch(bwa_ix_jobname,1,headpath(reference_path)) + bwaix_coms, debug)

    ## ------------------------------------------------- SBATCH COMMAND GENERATION  ---------------------------------------------------------------- ## 
    ## Initilizse list for sbatch
    sub_sbatchs = []
    ## ------------------------------------------------------- BWA INDEXING  ----------------------------------------------------------------------- ## 
    ## Submit the bwa index job
    if isbwaix(reference_path) or fileexists(bwaix_report):
        ## Print that we ahve detedted the indexed ref
        print('INFO: Detected bwa index, skipping indexing.')

        ## Set the job id 
        bwaix_jobid = ''

    else: ## Print to lien we are indexing 
        print('INFO: Submitting command to index reference with bwa.\nWARNING: Indexing may take a few hours to finish.')

        ## Format the in text
        intext = f'sbatch --dependency=singleton --partition={partition} {bwa_ix_jobname}'

        ## Append to sub list
        sub_sbatchs.append(intext)

        ## Print if we are in debug mode
        ifprint(intext,debug)

        ## Submit the bwa index 
        bwaix_jobid = fakejobid if (debug and runlocal) else submitsbatch(intext)

    ## Format the command file and initilize sample names
    command_files, samplenames, filteredbams = [(bwa_ix_jobname,'singleton',None,None,'bwaix',None,bwaix_report)], [], []

    ## Gather the fastqs 
    in_fastqs = getfastqs(fastqdir+'/')

    ## Assert we have fastq files
    assert len(in_fastqs) > 0, "ERROR: No fastq.gz files were detected!"

    ## iterate thru the pairs 
    for fastqix, (r1,r2) in enumerate(in_fastqs):
        ## Gather the sample name, exp mode, and fastp command file 
        sample_name =  getsamplename(r1)

        ## Append the sample names
        samplenames.append(sample_name)

        ## ------------------------------------------------------- FASTP SPLITTING ---------------------------------------------------------------- ## 
        ## Set the the fastp command file name 
        fastp_command_file =  f'{comsdir}/fastp.{sample_name}.sh'  

        fastp_coms, fastp_report = fastpeel(r1,r2,experi_mode,fastp_threads,fastp_splits)

        ## Write the command to file
        writetofile(fastp_command_file, sbatch(None,fastp_threads,the_cwd) + fastp_coms, debug)

        ## Append command to file
        command_files.append((fastp_command_file,'singleton',sample_name,experi_mode,'fastp','',fastp_report))

        ## Format the split fastq files
        fastq_splits = [(f'{splitsdir}/{formatn(i+1)}.{basename(r1)}', f'{splitsdir}/{formatn(i+1)}.{basename(r2)}') for i in range(fastp_splits) ]

        ## --------------------------------------------------------- ALIGNING -------------------------------------------------------------------- ## 
        ## Initlize list of bam files
        bam_files = []

        ## Format the bwa commands
        for (s1,s2) in fastq_splits:
            ## Call the bwa meme prep ftn, return the split bam file, the commands, and their report
            sbam, bwa_coms, bwa_report = prepbwamem(s1,s2,reference_path,bwa_threads,experi_mode)

            ## Format the bwa commands
            bwa_command_file =  f'{comsdir}/bwa.{basenobam(sbam)}.sh' 

            ## Write the bwa command to file 
            writetofile(bwa_command_file, sbatch(None,bwa_threads,the_cwd) + bwa_coms, debug)

            ## Append bwa command to list of commands 
            command_files.append((bwa_command_file,'afterok:',sample_name,experi_mode,'bwa',bwaix_jobid,bwa_report))

            ## ----------------------------------------------------- BAM SPLITTING ----------------------------------------------------------------- ## 
            ## Wrte command to split
            split_filename = f'{comsdir}/split.{basenobam(sbam)}.sh'

            ## Call the split ftn
            split_coms, split_report = splitcommand(sbam,mito,experi_mode,sam_threads,mapq)

            ## Write the split command
            writetofile(split_filename, sbatch(split_filename,sam_threads,the_cwd) + split_coms, debug)

            ## Apppend the seperating command
            command_files.append((split_filename,'afterok:',sample_name,experi_mode,'split','',split_report))

            ## Gather output bam name and txt names
            name_mapd_txt, name_plac_txt, name_mito_txt, name_unmp_txt, bedpe_txt = outnames(sbam,mito)    

            ## Remap the names to a bam file name 
            mapd_bam, plac_bam, mito_bam = txttobam(name_mapd_txt), txttobam(name_plac_txt), txttobam(name_mito_txt)
            
            ## Append the bam files
            bam_files.append((sbam, mapd_bam, plac_bam, mito_bam))

        ## --------------------------------------------------------- TXT and BEDPE MERGING ----------------------------------------------------------- ##
        ## Concat if we are in hic 
        if experi_mode == 'hic':
            for bamend in ['chrM.txt','placed.txt','bedpe']:
                ## Format command to merge the unmapped txt file
                concat_coms, concat_report = concatcom(f'{aligndir}/{sample_name}.{bamend}', f'{bamtmpdir}/*.{sample_name}.{bamend}')

                ## make concat file name
                concat_file = f'{comsdir}/concat.{sample_name}.{bamend}.sh'

                ## Write the concat command to file
                writetofile(concat_file, sbatch(concat_file,1,the_cwd) + concat_coms, debug)

                ## Append the concat command
                command_files.append((concat_file,'afterok:',sample_name,experi_mode,'merge','',concat_report))

        ## ----------------------------------------------------------- BAM MERGING ------------------------------------------------------------------ ##
        else: ## Format commands to merge on mapped, placed, and mitochondrial 
            for bamend in ['mapped','placed',mito]:
                ## Format the merge command 
                merge_coms, merge_report = mergebam(f'{sample_name}.{bamend}.bam', f'*.{sample_name}.{bamend}.bam', sam_threads)

                ## Format the name of the merge command 
                merge_file = f'{comsdir}/merge.{sample_name}.{bamend}.sh'

                ## Write the merge command to file 
                writetofile(merge_file, sbatch(merge_file,sam_threads,the_cwd) + merge_coms, debug)

                ## Append the merge command
                command_files.append((merge_file,'afterok:',sample_name,experi_mode,'merge','',merge_report))
        
        ## ------------------------------------------------------ MERGE UNMAPPED TXT ---------------------------------------------------------------- ##
        ## Format command to merge the unmapped txt file
        concat_coms, concat_report = concatcom(f'{aligndir}/{sample_name}.unmapped.txt', f'{bamtmpdir}/*.{sample_name}.unmapped.txt')

        ## make concat file name
        concat_file = f'{comsdir}/concat.{sample_name}.unmapped.sh'

        ## Write the concat command to file
        writetofile(concat_file, sbatch(concat_file,1,the_cwd) + concat_coms, debug)

        ## Append the concat command
        command_files.append((concat_file,'afterok:',sample_name,experi_mode,'merge','',concat_report))

        ## ---------------------------------------------------- MARKING DUPLICATES ----------------------------------------------------------------- ## 
        ## Format the merged bam file
        merged_map = f'{aligndir}/{sample_name}.mapped.bam'
        
        ## If we are skipping duplicates 
        if skipduplicates: ## Skip duplicate marking 
            ifprint('INFO: Skipping duplicate marking.',skipduplicates)
            
            ## Patch the merged to the marked bam 
            to_filt_bam = merged_map
            
        else: ## Otherwise call samblaster, make the file name 
            markdup_filename = f'{comsdir}/mark.{sample_name}.sh' 

            ## Call the mark duplicates ftn, returning the output name of the makred bam, the commands for makring, and the assoicated report 
            to_filt_bam, mark_coms, mark_report = markduplicates(merged_map,sam_threads)

            ## Write to file the mark dups commands
            writetofile(markdup_filename, sbatch(markdup_filename,sam_threads,the_cwd) + mark_coms, debug)

            ## Append the marking duplicates command 
            command_files.append((markdup_filename,'afterok:',sample_name,experi_mode,'mark','',mark_report))

        ## ---------------------------------------------------- FILTER BAM FILE ------------------------------------------------------------------ ## 
        ## Format filtering command
        filter_filename = f'{comsdir}/filter.{sample_name}.sh'

        ## Call the filter command 
        filt_bam, filt_coms, filt_report = filterbam(to_filt_bam,mapq,sam_threads)

        ## Write to file the filt commands and sbatch
        writetofile(filter_filename, sbatch(filter_filename,sam_threads,the_cwd) + filt_coms, debug)

        ## Append the filtering command
        command_files.append((filter_filename,'afterok:',sample_name,experi_mode,'filter','',filt_report))

        ## Append the output filtered bam
        filteredbams.append(filt_bam)
    ## ------------------------------------------------------------------------------------------------------------------------------------------ ##
    ## <- This indent marks all the following commands run post aligment and filtering of individual samples 
    ## ------------------------------------------------- COUNTING BAM FILE ---------------------------------------------------------------------- ## 
    ## Format counting command, the count report, and the command
    counting_filename, count_report = f'{comsdir}/countbams.{sample_name}.sh', reportname(f'{aligndir}/bam.counts.bam','count')

    ## Format the save file name
    save_dist_name = f'-S ./{diagdir}/{run_name}.bam.fragment.dist.png\n'

    ## Add fragment diagnostic plots 
    if experi_mode == 'atac':
        ## Format the fragment histogram distribution for atac seq samples 
        frag_calc_commands = f'{scriptsdir}/pymacs2.py -r {run_name}\n{scriptsdir}/fragmentdist.py -b ./{aligndir}/*.primary.*.bam ' + save_dist_name
    ## Add fragment diagnostic plots for chip mode 
    elif experi_mode == 'chip':
        ## Format the contorl samples
        joined_controls = ' '.join(chip_control)
        ## Format the fragment commands
        frag_calc_commands = f'{scriptsdir}/pymacs2.py -r {run_name}\n{scriptsdir}/fragmentdist.py -b ./{aligndir}/*.primary.*.bam {joined_controls} ' + save_dist_name
    else: ## Otherwise pass a new line character 
        frag_calc_commands = '\n'

    ## List the count commands 
    count_commands = [f'{scriptsdir}/countbams.py {run_name} {experi_mode}\n', frag_calc_commands, f'echo Finished counting bam files in {aligndir} dir. >> {count_report}\n']

    ## Wriet the coutn command to file
    writetofile(counting_filename,sbatch(counting_filename,1,the_cwd) + count_commands, debug)

    ## Append the counting command
    command_files.append((counting_filename,'afterok:',sample_name,experi_mode,'count','',count_report))

    ## ---------------------------------------------------- PEAK CALLING MACS2 ----------------------------------------------------------------- ## 
    ## If we are running analysis on chip or atac experiments 
    if experi_mode in ['atac','chip']:
        ## Gather the genome size 
        gsize = genomesize(g_size,reference_path,mito)

        ## Format the macs2 call report name
        macs2_report, macs2_filename = reportname(run_name,'macs2'), f'{comsdir}/macs2.{run_name}.sh'
    
        ## Format the command to macs2
        macs2_commands = peakattack(filteredbams,run_name,macs2_report,incontrols=chip_control,gsize=gsize,broad=broadpeak)

        ## Write the macs2 commands to file
        writetofile(macs2_filename, sbatch(macs2_filename,1,the_cwd) + macs2_commands, debug)

        ## Append the macs2 command 
        command_files.append((macs2_filename,'afterok:',sample_name,experi_mode,'macs2','',macs2_report))
    else:
        pass 
    ## ------------------------------------------------- SUBMITTING TIME-STAMP ----------------------------------------------------------------- ## 
    ## Set the timesampe file name 
    timestamp_file, timestampsh, timestamp_report = f'{diagdir}/{run_name}.timestamp.{stamp}.txt', f'{comsdir}/time.stamp.sh', reportname(run_name,f'timestamp.{stamp}')

    ## Formath time stamp and echo commands 
    times_commands = [f'{scriptsdir}/endstamp.py {timestamp_file} {stamp}\n', f'echo Finished SLURPY run of sample: {run_name}. >> {timestamp_report}\n']

    ## Format the command file name and write to sbatch, we will always ask the timestamp to run even in debug mode 
    writetofile(timestampsh, sbatch(timestampsh,1,the_cwd) + times_commands, False)

    ## Append the timestamp command to file
    command_files.append((timestampsh,'afterok:',run_name,experi_mode,'timestamp','',timestamp_report))

    ## ----------------------------------------------------- CLEANING UP FILES ------------------------------------------------------------------- ## 
    if ifclean or (rerun == 'clean'): 
        ## Format command to remove uneedeed files at end of run 
        remove_sh, remove_report = f'{comsdir}/cleanup.sh', reportname(run_name,'clean')

        ## Format the command to clean up          
        # ./{aligndir}/*primary*\n'] # <- remove primary alignments. We are not doing this anymore. 
        writetofile(remove_sh, sbatch(remove_sh,1,the_cwd) + [f'rm -rf {bamtmpdir} {splitsdir}'], debug)

        ## Append the clean up command to file
        command_files.append((remove_sh,'afterok:',run_name,experi_mode,'clean','',remove_report))
    else: ## Otherwise do nothing
        pass 
    ## --------------------------------------------------------- RESTARTING  -------------------------------------------------------------------- ## 
    ## Call the command dataframe
    command_files, was_hard_reset = commandcontrol(command_files,hardreset,pipeline_steps,rerun,bwaix_jobid)

    ## ---------------------------------------------------- COMMAND SUBMISSION ------------------------------------------------------------------ ## 
    ## ------------------------------------------------ SUBMITTING FASTP SPLITTING -------------------------------------------------------------- ## 
    ## Iterate thru the fastp opperations 
    for rix,row in command_files[(command_files.Operation=='fastp') & (command_files.Torun==0)].iterrows():

        ## Format the jobparts
        jobparts = jobname(['fastp',rix,partition,stamp])

        ## Format the sbatch   
        intext = f'sbatch --dependency=singleton {jobparts} --partition={partition} {row.Jobfile}'

        ## Append to sub list
        sub_sbatchs.append(intext)

        ## Append the fastp job id 
        command_files.loc[rix,'JobID'] = fakejobid + rix if (debug and runlocal) else submitsbatch(intext)

        ## Print if we are in debug mode
        ifprint(intext,debug)

    ## --------------------------------------------------- SUBMITTING ALIGNMENTS --------------------------------------------------------------- ## 
    ## Gather the fastp job ids
    fastpids = formatids(command_files,['fastp','bwaix']) 

    ## Gather the bwa index 
    bwa_indexs = command_files[(command_files.Torun==0) & (command_files.Operation=='bwa')].index.tolist()

    ## Group the bwa commands 
    grouped_bwa = [bwa_indexs[i::bwa_runs] for i in range(bwa_runs)]

    ## Iteraet thru the groups of the bwa comamnds, then the groups
    for gi, gbwa in enumerate(grouped_bwa):
        for i, bwai in enumerate(gbwa):
            ## Set the row fromt he command files 
            row = command_files.loc[bwai,:]

            ## Format the job parts
            jobparts = jobname(['bwa',gi,partition,stamp])

            ## If this is the first in the bwa submission, add the fastp dependency, otherwise make the bwa a named singleton 
            intext = f'sbatch{formatdepends(fastpids)}{jobparts} --partition={partition} {row.Jobfile}' if (i==0) else f'sbatch --dependency=singleton {jobparts} --partition={partition} {row.Jobfile}'

            ## Append to sub list
            sub_sbatchs.append(intext)

            ## Submit the job and append the job id to the dataframe 
            command_files.loc[bwai,'JobID'] = fakejobid + bwai if (debug and runlocal) else submitsbatch(intext)

            ## Print if we are in debug mode
            ifprint(intext,debug)

    ## ---------------------------------------------- SUBMITTING SPLITER COMMANDS ----------------------------------------------------------- ## 
    ## Submit the splitting commands and append the submitted jobs
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'split','bwa',stamp,partition,bylast=True)

    ## ---------------------------------------------- SUBMITTING MERGE COMMANDS ------------------------------------------------------------- ## 
    ## Submit the merge commands and append to sub 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'merge','split',stamp,partition)

    ## ------------------------------------------ SUBMITTING MARK DUPLICATES COMMAND -------------------------------------------------------- ## 
    ## Submit the mark duplicate commands and append to sub 
    sub_mark = [] if skipduplicates else submitdependency(command_files,'mark','merge',stamp,partition)
    ## Add the sbatchs 
    sub_sbatchs = sub_sbatchs + sub_mark 

    ## ---------------------------------------------- SUBMITTING FILTER COMMANDS ------------------------------------------------------------ ## 
    ## Set the above step 
    above_step = 'merge' if skipduplicates else 'mark'
    ## Sumit the filtering commands and submit to file
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'filter',above_step,stamp,partition)
    
    ## ----------------------------------------------- SUBMITTING MACS2 COMMANDS ------------------------------------------------------------ ## 
    ## Submit call to macs2 if in atac or chip mode
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'macs2','filter',stamp,partition,group='Experiment')

    ## ----------------------------------------------- SUBMITTING COUNT COMMANDS ------------------------------------------------------------ ## 
    ## Submit the count command 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'count','filter',stamp,partition,group='Experiment')

    ## ------------------------------------------------ SUBMITTING TIME COMMANDS ------------------------------------------------------------ ## 
    ## Set the above step
    above_step = 'macs2' if experi_mode in ['atac','chip'] else 'filter'
    ## Submit time stamp 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'timestamp',[above_step,'count'],stamp,partition,group='Experiment')

    ## ---------------------------------------------------- CLEAN UP COMMANDS --------------------------------------------------------------- ## 
    ## Submit the clean up command if the flag was passed
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'clean','timestamp',stamp,partition,group='Experiment')

    ## ----------------------------------------------- SAVING OUT FILE & COMMANDS ----------------------------------------------------------- ## 
    ## Write the sbatch commands
    writetofile(f'{debugdir}/sbatch.log.txt',[l+'\n' for l in sub_sbatchs], False)

    ## Calc number of jobs submitted
    njobstosub, njobssubbed = command_files[(command_files.Torun==0)].shape[0], len(sub_sbatchs)

    ## Check our work
    ifprint(f'WARNING: The number of expected jobs to run ({njobstosub}) and number of jobs submitted {njobssubbed} does not match!' ,not (njobstosub == njobssubbed))

    ## Print the number of commands being submitted
    print(f'INFO: A total number of {njobssubbed} jobs were submitted with this run.')

    ## If zero jobs were submitted
    ifprint('WARNING: Zero jobs were submitted; If this was unexpected try running slurpy again, including --restart flag.',(njobssubbed==0))

    ## Save out the command files 
    command_files.to_csv(f'{debugdir}/command.file.csv', index=False)
## End of file 