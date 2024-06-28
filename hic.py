#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      hic.py; A Hi-C Processing Pipeline 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## List command for canceling all
"""
squeue -u croth | grep 'croth' | awk '{print $1}' | xargs -n 1 scancel
squeue -u croth | grep 'croth' | grep tb | awk '{print $1}' | xargs -n 1 scancel
squeue -u croth | grep 'croth' | grep gpu | grep "(DependencyNeverSatisfied)" | awk '{print $1}' | xargs -n 1 scancel

""" 
## List slurpy command for VERO analysis
"""
./SLURPY/hic.py -r /panfs/biopan04/epier/Cryo_ATAC/GreenMVA_Ref/Chlorocebus_sabeus_mva.fasta -M NC_008066.1 -G ../Chlorocebus_sabeus_mva.genome.sizes.autosome.filtered.bed --merge 

"""

##      MODULE LOADING and VARIABLE SETTING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in ftns and variables from defaults
from defaults import *
## Load in ftns from other libraries
from pysamtools import checksam, writetofile, chromdf
## Bring in bwa mem ftn for hic
from pybwatools import bwamem_hic
 
## Set the ftn descritption
hiclite_descr = "Processing and analysis pipeline for paired-end sequencing data from Hi-C experiments."
## Set the hardset hic mode 
experi_mode = 'hic'
## Define help messages
R_help = "Step within the pipeline to re-run from. Options for Hi-C analysis include: %s"%(', '.join(hic_pipeline))
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Ftn for formating fastp command to filter and split reads
def hicpeel(r1:str, r2:str, w:int, s:int, z=4, options=fastp_opts, script='fastp') -> tuple:
    """Formats calls to fastp given input fastq files."""
    ## Format the splits
    split1, split2   = f'{splitsdir}/{basename(r1)}', f'{splitsdir}/{basename(r2)}'  
    ## Set the report name
    report = reportname(getsamplename(r1),script)
    ## Gather the fastp header, json file, and html 
    fastp_header,fastp_json,fastp_html = fastcut(r1,r2,split1,split2,report,1)
    ## Call fastp again to split the files
    slice = fastp_header + f' -s {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'
    ## Call the fast dry command and format the remove command 
    dry = fastdry(r1,r2,report)
    ## Format and return the command based on the experiment type (ie if it is hi-c or not)
    return [slice,dry,f'mv {fastp_html} {diagdir}/\n',f'mv {fastp_json} {diagdir}/\n'], report 

## Ftn for formating fastq files for bwa mem in hic mode
def prepbwamem(r1:str, r2:str, refix:str, t:int) -> tuple:
    """Formats a bwa mem command given the paired reads, reference index, thread count, and mode."""
    ## Format the sam file name and report name
    sam, report = f'{bamtmpdir}/{getsamplename(r1)}.sam', f'{debugdir}/{getsamplename(r1)}.bwa.log.txt'
    ## Format the bwa mem command 
    makebams = bwamem_hic(r1,r2,refix,sam,report,threads=t)
    ## Return the sam, make bam command and reprot 
    return sam, makebams, report  

## Ftn for splitting hic files by chromosmoe
def hicspliter(inpath:str ,outpath:str, chrom:str, chromosomes:list, dedup:bool, script='split') -> tuple:
    ## Format the report name
    report = reportname(inpath,script)
    ## Set the remove duplicates argument 
    remove_du = ' --skip-dedup' if dedup else ''
    ## Set the pathtochroms
    chromosomes = ' '.join(chromosomes)
    ## List the commands 
    command_lines = [f'{scriptsdir}/splithic.py -i {inpath} -o {outpath} -c {chrom} -g {chromosomes}{remove_du}\n', 
                     f'echo Finished duplicate filtering and spliting Hi-C contacts by chromosome for {chrom} from {inpath}. >> {report}\n']
    ## Return the command lines and report  
    return command_lines, report 

## Ftn for sorting hic file
def hicsorter(inpaths:list, outname:str, script='sort') -> tuple:
    ## Format the report name
    report = reportname(inpaths[0],script)
    ## Format and return the call to the hic sort command and report echo 
    command_lines = [f'{scriptsdir}/sorthic.py -i ' +' '.join(inpaths) + f' -o {outname}\n', 
                     f'echo Finished sorting Hi-C contacts by chromosome from {len(inpaths)} files. >> {report}\n']
    ## Return the commands and the report 
    return command_lines, report 

## Ftn for calling juicer's pre command
def juicerpre(intxt:str, outhic:str, Xmemory:int, jarfile:str, threadcount:int, bins:list, genomepath:str, script='juicerpre') -> tuple:
    """
    java -Xmx49152m -Xms49152m -jar $jarpath pre -j 5 -r 500000,250000,200000,150000,100000,50000,25000,10000,5000,1000 $1 ${1}.hic $2
    """
    ## Format the report name 
    report = reportname(outhic+'.bam',script)
    ## Set the java command for the passed juicer jar file 
    prestr = ['java -Xmx%sm -Xms%sm -jar %s pre -j %s -r %s %s %s %s\n'%(Xmemory,Xmemory,jarfile,threadcount,','.join(map(str,bins)),intxt,outhic,genomepath),
              f'echo Finished formating Hi-C contacts into a .hic file on path: {outhic} >> {report}\n']

    ## Return the pre and report
    return prestr, report

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is envoked by name 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = hiclite_descr)

    ## Add the required argument
    parser.add_argument("-r", "--refix",          dest="r", type=str,  required=True, help = r_help,  metavar = './path/to/reference.bwaix'                    ) 
    
    ## Add the default arguments
    parser.add_argument("-F", "--fastp-splits",   dest="F", type=int,  required=False, help = F_help, metavar = splitsize,              default = splitsize    )
    parser.add_argument("-B", "--parallel-bwa",   dest="B", type=int,  required=False, help = B_help, metavar = parallelbwa,            default = parallelbwa  )
    parser.add_argument("-P", "--partition",      dest="P", type=str,  required=False, help = P_help, metavar = part,                   default = part         ) 
    parser.add_argument("-M", "--mtDNA",          dest="M", type=str,  required=False, help = M_help, metavar = mito,                   default = mito         )
    parser.add_argument("-X", "--exclude",        dest="X", nargs='+', required=False, help = X_help, metavar = 'chrX, chrY ...',       default = []           )
    parser.add_argument("-Q", "--map-threshold",  dest="Q", type=int,  required=False, help = Q_help, metavar = map_q_thres,            default = map_q_thres  )
    parser.add_argument("-R", "--rerun-from",     dest="R", type=str,  required=False, help = R_help, metavar = 'step',                 default = None         )
    parser.add_argument("-q", "--fastq",          dest="q", type=str,  required=False, help = q_help, metavar = '.fastq.gz',            default = fends        )

    ## Set number of threads 
    parser.add_argument("-f", "--fastp-threads",  dest="f", type=int,  required=False, help = f_help, metavar = fastpthreads,           default = fastpthreads )
    parser.add_argument("-b", "--bwa-threads",    dest="b", type=int,  required=False, help = b_help, metavar = bwathreads,             default = bwathreads   )
    parser.add_argument("-t", "--dask-threads",   dest="t", type=int,  required=False, help = t_help, metavar = daskthreads,            default = daskthreads  )

    ## Set values for Hi-C analysis 
    parser.add_argument("-n", "--run-name",       dest="n", type=str,  required=False, help = n_help, metavar = 'name',                 default = None         )
    parser.add_argument("-E", "--error-distance", dest="E", type=int,  required=False, help = E_help, metavar = 'bp',                   default = error_dist   )
    parser.add_argument("-C", "--self-circle",    dest="C", type=int,  required=False, help = C_help, metavar = 'bp',                   default = circle_dist  )
    parser.add_argument("-L", "--library",        dest="L", type=str,  required=False, help = L_help, metavar = 'MboI',                 default = lib_default  )
    parser.add_argument("-D", "--mindist",        dest="D", type=int,  required=False, help = D_help, metavar = 'n',                    default = set_distance )
    parser.add_argument("-Z", "--chunksize",      dest="Z", type=int,  required=False, help = Z_help, metavar = 'n',                    default = chunks       )
    parser.add_argument("-G", "--genomelist",     dest="G", type=str,  required=False, help = G_help, metavar = './path/to/list.tsv',   default = None         )
    parser.add_argument("-J", "--jar-path",       dest="J", type=str,  required=False, help = J_help, metavar = './path/to/juicer.jar', default = None         )
    parser.add_argument("-x", "--Xmemory",        dest="x", type=int,  required=False, help = x_help, metavar = xmemory,                default = xmemory      )
    parser.add_argument("-S", "--bin-sizes",      dest="S", nargs='+', required=False, help = S_help, metavar = '25000, 10000, ...',    default = binsizes     )

    ## Set boolean flags 
    parser.add_argument("--restart",              dest="start",  help = restart_help,  action = 'store_true')
    parser.add_argument("--debug",                dest="debug",  help = debug_help,    action = 'store_true')
    parser.add_argument("--skip-dedup",           dest="mark",   help = mark_help,     action = 'store_true')
    parser.add_argument("--clean",                dest="clean",  help = clean_help,    action = 'store_true')
    parser.add_argument("--merge",                dest="merge",  help = merge_help,    action = 'store_true')


    ## Set the paresed values as inputs
    inputs = parser.parse_args()
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    

    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set required variables
    reference_path  = inputs.r       ##     Set path to the reference genome
                                     ##   
    ## Set default vairables         ##
    fastp_splits    = inputs.F       ##     Number of splits in fastp 
    bwa_runs        = inputs.B       ##     Set the number of parallel runs of bwa 
    partition       = inputs.P       ##     Set the partition 
    mito            = inputs.M       ##     Set the mito contig name 
    excludes        = inputs.X       ##     List of chromosomes to exclude from analysis 
    mapq            = inputs.Q       ##     Set the mapping quality threshold 
    rerun           = inputs.R       ##     Setp to rerun pipeline from 
    fend            = inputs.q       ##     End of the input fastq files 
                                     ##    
    ## Set threads                   ##   
    fastp_threads   = inputs.f       ##     Number of fastp threads
    bwa_threads     = inputs.b       ##     Number of threads in bwa alignments
    dask_threads    = inputs.t       ##     Number of threads in dask
                                     ##       
    ## Set variables for Hi-C        ## 
    run_name        = inputs.n       ##     The name of the samples 
    error_distance  = inputs.E       ##     Distance in bp to check hi-C contacts for errors 
    circle_distance = inputs.C       ##     Distacne in bp to check for self-circle
    enzyme          = inputs.L       ##     Restriction-enzyme used in Hi-C prep 
    set_distance    = inputs.D       ##     Overall minimum intra-chromosomal contact distance allowed 
    chunk_size      = inputs.Z       ##     Chunk size (row number) to load in with pandas 
    pathtochrom     = inputs.G       ##     Path to list of chromosomes to use 
    jarpath         = inputs.J       ##     Path to juicer jar file 
    xmemory         = inputs.x       ##     Amount of memory passed to juicer pre command 
    binsizes        = inputs.S       ##     Bins / resolutions used in hi-c analysis 
                                     ##
    ## Set boolean vars              ##
    hardreset       = inputs.start   ##     Resetart the slurpy run, removing previous
    debug           = inputs.debug   ##     Run in debug mode 
    skipduplicates  = inputs.mark    ##     Boolean to mark duplicates 
    ifclean         = inputs.clean   ##     Flag to run clean up script 
    postmerging     = inputs.merge   ##     Forces premerge of outputs 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      ROTH SETTINGS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Reset reference with presets if we are running as Cullen Roth (These can be edited if you are not me :-))
    reference_path = t2t_refpath if reference_path in ['T2T','t2t','T2t','t2T'] else reference_path
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      INITILIZATION 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Check that the fastq path exists
    assert pathexists(fastqdir), "ERROR: Unable to detect a fastqs directory!"
    ## Check the versions of samtools, the user email is an email and the experiment mode is one we know
    assert checksam(), not_sam_err 
    ## If needed reset that the fastp threads and splits such that they are a multiple of each
    fastp_splits,fastp_threads = checkfastp(fastp_splits,fastp_threads)
    ## Reformat clean boolean if clean was passed from restart
    ifclean = patchclean(rerun,ifclean)
    ## Set the hic file ends
    hicfileends = [mito] + hicfileends_tmp
    ## Initilizse list for sbatch
    sub_sbatchs = []

    ## Set the lib error
    lib_error = "ERROR: The passed enzyme used in Hi-C library prep was not recognized."
    ## Check library if it was pased
    if enzyme:
        assert enzyme.lower() in ['mboi', 'dpnii','sau3ai','hindiii','arima','none'], lib_error

    ## Reset the library enzyme variable if the none option was passed
    #enzyme = None if enzyme.lower() == 'none' else enzyme ## We don't want to do this 

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    ##      CONFIRM THE HARD RESTART 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## If there is a hard reset passed 
    confirmreset(grouped_dirs) if hardreset else None
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      MORE MODULE LOADING (and time stamp setting)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in the list dir ftn 
    from os import getcwd as gcwd
    ## Load in time mod
    import time 
    ## Set the unique time stamp
    stamp = time.time()





    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      DIRECTORY MAKING & TIME STAMP SUBMISSION 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Let the user know we are making directories 
    print(directormaking)
    ## Get the current working dir
    the_cwd = gcwd()
    ## If the run name is none
    run_name = setrunname(run_name,the_cwd)
    #print(run_name)
    ## Make the directories 
    makedirectories(grouped_dirs)
    ## Gather and remove the old command files, error, and output logs if they are there 
    [remove(f) for f in sortglob(f'{comsdir}/*.sh') + sortglob(f'{debugdir}/*.err') + sortglob(f'{debugdir}/*.out')]
    ## Generate the bwa indx command
    bwa_ix_jobname = f'{comsdir}/bwaindex.sh'
    ## Call the bwa command and its report 
    bwaix_coms, bwaix_report = bwaindex(reference_path)
    ## Write the bwa index command and sbatch to file 
    writetofile(bwa_ix_jobname, sbatch(bwa_ix_jobname,1,headpath(reference_path)) + bwaix_coms, debug)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    ##      CHROMOSOME GATHERING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Inform the user we are gathering chromosomes
    print(chromgathering)
    ## Set fai path
    fai_path = reference_path + '.fai'
    ## Set new path to chrom
    newpathtochrom = aligndir + '/' + reference_path.split('/')[-1].split('.fa')[0] + '.txt'
    ## If a path of chrom file was passed and exists 
    if pathtochrom and fileexists(pathtochrom): 
        ## Set the list of chromosomes 
        chrlist = readtable(pathtochrom)[0].values
    elif fileexists(fai_path):
        ## Patch path
        pathtochrom = newpathtochrom
        ## Load in the fai file fromt he reference 
        chrbed = readtable(fai_path)
        ## SAve out the bed file 
        chrbed[[0,1]].to_csv(pathtochrom,sep=' ',index=False,header=False)
        ## Set the list of chromosomes 
        chrlist = chrbed[0].values
    ## Make a dataframe of the files 
    else: ## Bring in the list of chromosomes from the fasta file
        ## Patch chrom 
        pathtochrom = newpathtochrom
        ## IF, the path to chrom form the reference file path already exists 
        if fileexists(pathtochrom):
            chrbed = readtable(pathtochrom)
        else:
            ## Gather tupes of chromosome ids and lengths 
            chrbed = chromdf(reference_path)
            ## Save out the chrbed as an .txt file for next time or further analysis 
            chrbed.to_csv(pathtochrom,sep=' ',index=False,header=False)
        ## Set the list of chromosomes 
        chrlist = chrbed[0].values

    ## Expand exlcude list to include mitochondria contic 
    excludes.append(mito)
        
    ## Remove mito
    chrlist = [c for c in chrlist if (c not in excludes)]
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      BWA INDEXING  
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
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
        bwaix_jobid = fakejobid if debug else submitsbatch(intext)

    ## Format the command file and initilize sample names
    command_files, samplenames, filteredbams = [(bwa_ix_jobname,'singleton',None,None,'bwaix',None,bwaix_report)], [], []
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      FASTQ GATHERING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Inform user we are formating jobs
    print(formatingfastq)
    ## Gather the fastqs 
    in_fastqs = getfastqs(fastqdir+'/*.gz')
    ## Assert we have fastq files
    assert len(in_fastqs) > 0, "ERROR: No fastq.gz files were detected!"
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SAMPLE NAME MAKING AND HANDELING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## iterate thru the pairs 
    for fastqix, (r1,r2) in enumerate(in_fastqs):
        ## Gather the sample name, exp mode, and fastp command file 
        sample_name = getsamplename(r1) 
        ## Append the sample names
        samplenames.append(sample_name)
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      FASTP SPLITTING 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ##  Set the the fastp command file name 
        fastp_command_file =  f'{comsdir}/fastp.{sample_name}.sh'  
        ##  Gather the fastp command and report 
        fastp_coms, fastp_report = hicpeel(r1,r2,fastp_threads,fastp_splits)
        ##  Write the command to file
        writetofile(fastp_command_file, sbatch(None,fastp_threads,the_cwd) + fastp_coms, debug)
        ##  Append command to file
        command_files.append((fastp_command_file,'singleton',sample_name,experi_mode,'fastp','',fastp_report))
        ##  Format the split fastq files
        fastq_splits = [(f'{splitsdir}/{formatn(i+1)}.{basename(r1)}', f'{splitsdir}/{formatn(i+1)}.{basename(r2)}') for i in range(fastp_splits)]
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      PARALLIZATION of ALIGNMENT and Hi-C Filtering 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Format the bwa commands
        for (s1,s2) in fastq_splits:
            ## Call the bwa meme prep ftn, return the split bam file, the commands, and their report
            sbam, bwa_coms, bwa_report = prepbwamem(s1,s2,reference_path,bwa_threads)
            ## Format the bwa commands
            bwa_command_file =  f'{comsdir}/bwa.{basenobam(sbam)}.sh' 
            ## Write the bwa command to file 
            writetofile(bwa_command_file, sbatch(None,bwa_threads,the_cwd) + bwa_coms, debug)
            ## Append bwa command to list of commands 
            command_files.append((bwa_command_file,'afterok:',sample_name,experi_mode,'bwa',bwaix_jobid,bwa_report))
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##


            ##      PRE FILTERING (post bwa) 
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## Format prefilter report name, commands, and file name 
            pre_filt_repo = reportname(sbam,'pre')                                                  ##      Generate report file name
            pre_filt_coms = [f'{scriptsdir}/prefilter.py -s {sbam} -M {mito} -Z {chunk_size}\n',    ##      Format the commands for pre filtering
                             f'echo Finished prefiltering on {sbam} >> {pre_filt_repo}']            ##  
            pre_filt_file = f'{comsdir}/prefilt.{basenobam(sbam)}.sh'                               ##      Set teh pre filt bash file name 
            ## Write the above to file
            writetofile(pre_filt_file, sbatch(pre_filt_file,2,the_cwd)+pre_filt_coms, debug)
            ## Append to command file
            command_files.append((pre_filt_file,'afterok:',sample_name,experi_mode,'pre','',pre_filt_repo))
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##


            ##      POST FILTERING 
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## Format the name of the genomic mappings filtered from above
            tsam = sbam.split('.sam')[0] + '.genomic.txt'
            ## Format postfilter report name, commands, and file name 
            post_filt_repo = reportname(sbam,'post')
            post_filt_coms = [f'{scriptsdir}/postfilter.py -t {tsam} -r {reference_path} -Q {mapq} -D {set_distance} -Z {chunk_size}\n', 
                              f'echo Finished postfiltering on {tsam} >> {post_filt_repo}']
            post_filt_file = f'{comsdir}/postfilt.{basenobam(tsam)}.sh'
            ## Write the above to file
            writetofile(post_filt_file, sbatch(post_filt_file,2,the_cwd)+post_filt_coms, debug)
            ## Append to command file
            command_files.append((post_filt_file,'afterok:',sample_name,experi_mode,'post','',post_filt_repo))
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##


            ##      Hi-C FILTERING 
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## Format the expected text file name for post filtering 
            fsam = tsam.split('.txt')[0] + '.filtered.txt'
            ## Format hic-filter report name, commands, and file name 
            hic_filt_repo = reportname(sbam,'filter')
            hic_filt_coms = [f'{scriptsdir}/filterhic.py -t {fsam} -r {reference_path} -l {enzyme} -C {circle_distance} -E {error_distance} -Z {chunk_size}\n',
                             f'echo Finished postfiltering on {fsam} >> {hic_filt_repo}']
            hic_filt_file = f'{comsdir}/hicfilt.{basenobam(fsam)}.sh'
            ## Write the above to file
            writetofile(hic_filt_file, sbatch(hic_filt_file,2,the_cwd)+hic_filt_coms, debug)
            ## Append to command file
            command_files.append((hic_filt_file,'afterok:',sample_name,experi_mode,'filter','',hic_filt_repo))
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      FILE and TXT MERGING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Iterate over the sample names
    for si,sname in enumerate(samplenames):
        ## Set the smaple name 
        sample_name = run_name if postmerging else sname 
        ## If the sample is the first and we are merging all smaples 
        if (si > 0) and postmerging:
            break 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

        
        ##      TEXT (HIC) MERGING
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Iterate over the end file names used in hi-C analysis, except the last tohic files.
        for i,bamend in enumerate(hicfileends):
            ## Format the infilename, we need this last filename in the loop for the hic
            newcatfile =  f'{aligndir}/{sample_name}.{bamend}.txt'
            ## Format the outfile 
            wildcard_filename = f'{bamtmpdir}/*.{bamend}.txt' if postmerging else f'{bamtmpdir}/*.{sample_name}*.{bamend}.txt'
            ## Format command to merge the unmapped txt file
            concat_coms, concat_report = pandacat(wildcard_filename,newcatfile,rmheader=(bamend == 'tohic')) 
            ## make concat file name
            concat_file = f'{comsdir}/concat.{sample_name}.{bamend}.sh'
            ## Write the concat command to file
            writetofile(concat_file, sbatch(concat_file,1,the_cwd) + concat_coms, debug)
            ## Append the concat command
            command_files.append((concat_file,'afterok:',sample_name,experi_mode,'concat','',concat_report))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        

        ##      SPLITING and MARKING DUPLICATES in Hi-C FILES
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Initilizse list of chromosome calls 
        chrom_cats = [] 
        ## Iterate thru each chromosome
        for chrom in chrlist:
            ## Set the chrom cat file
            chromcatfile =  f'{bamtmpdir}/{sample_name}.{bamend}.{chrom}.txt' 
            ## append file name
            chrom_cats.append(chromcatfile)
            ## Sort, duplicate remove, and saveout the the concatonated hic file 
            split_coms, split_report = hicspliter(newcatfile,chromcatfile,chrom,chrlist,skipduplicates)
            ## make concat file name
            split_file = f'{comsdir}/split.{sample_name}.{bamend}.{chrom}.sh'
            ## Write the concat command to file
            writetofile(split_file, sbatch(split_file,2,the_cwd) + split_coms, debug)
            ## Append the concat command
            command_files.append((split_file,'afterok:',sample_name,experi_mode,'split','',split_report))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      SORT the Hi-C Contacts 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Set the output file 
        outtxtpath = newcatfile.split('.txt')[0] + ('.sorted.txt' if skipduplicates else '.nodups.sorted.txt')
        ## Print if we are skipping duplicates
        ifprint('INFO: Skipping duplicate marking.',skipduplicates)

        ## Sort and saveout the the concatonated hic file 
        sort_coms, sort_report = hicsorter(chrom_cats,outtxtpath)
        ## make concat file name
        sort_file = f'{comsdir}/sort.{sample_name}.{bamend}.sh'
        ## Write the concat command to file
        writetofile(sort_file, sbatch(sort_file,dask_threads,the_cwd) + sort_coms, debug)
        ## Append the concat command
        command_files.append((sort_file,'afterok:',sample_name,experi_mode,'sort','',sort_report))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      Make .hic file(s) with Juicer's pre command 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Set if jar was passed
        jargiven = jarpath and fileexists(jarpath)
        ## IF we have a valid jar path
        if jargiven:
            ## Set the output hic path
            outhicpath = outtxtpath.split('.txt')[0] + '.hic'
            ## Sort and saveout the the concatonated hic file 
            jpre_coms, jpre_report = juicerpre(outtxtpath,outhicpath,xmemory,jarpath,fastp_threads,binsizes,pathtochrom,fastp_threads)
            ## make concat file name
            jpre_file = f'{comsdir}/juicerpre.{sample_name}.{bamend}.sh'
            ## Write the concat command to file
            writetofile(jpre_file, sbatch(jpre_file,fastp_threads,the_cwd) + jpre_coms, debug)
            ## Append the concat command
            command_files.append((jpre_file,'afterok:',sample_name,experi_mode,'juicerpre','',jpre_report))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## <- This indent marks all the following commands run post aligment and filtering of individual samples 
    

    ##      COUNTING TXT FILES      
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Format counting command, the count report, and the command
    counting_filename = f'{comsdir}/count.hic.{run_name}.sh'         ##     Name the bash command file 
    count_report      = reportname(f'{aligndir}/hic.bam','count')    ##     Set the report name 
    ## List the count commands 
    count_commands = [f'{scriptsdir}/counthic.py {run_name}\n', f'echo Finished counting bam files in {aligndir} dir. >> {count_report}\n']
    ## write the count command to file
    writetofile(counting_filename,sbatch(counting_filename,1,the_cwd) + count_commands, debug)
    ## append the counting command
    command_files.append((counting_filename,'afterok:',run_name,experi_mode,'count','',count_report))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SUBMITTING TIME-STAMP   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set the timesampe assocaited file names 
    timestamp_file   = f'{diagdir}/{run_name}.timestamp.{stamp}.txt' ##     Name of the output file 
    timestampsh      = f'{comsdir}/time.stamp.sh'                    ##     Name of the .sh bash file 
    timestamp_report = reportname(run_name,f'timestamp.{stamp}')     ##     Name of the log to report to 
    ## Formath time stamp and echo commands 
    times_commands = [f'{scriptsdir}/endstamp.py {timestamp_file} {stamp}\n', f'echo Finished SLURPY run of sample: {run_name}. >> {timestamp_report}\n']
    ## Format the command file name and write to sbatch, we will always ask the timestamp to run even in debug mode 
    writetofile(timestampsh, sbatch(timestampsh,1,the_cwd) + times_commands, False)
    ## Append the timestamp command to file
    command_files.append((timestampsh,'afterok:',run_name,experi_mode,'timestamp','',timestamp_report))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CLEANING UP FILES   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## If the clean boolean or rerun vars were set 
    if ifclean or (rerun == 'clean'): 
        ## Format command to remove uneedeed files 
        remove_sh     = f'{comsdir}/cleanup.sh'        ##   Set the bash file name 
        remove_report = reportname(run_name,'clean')   ##   Set the report 
        ## Format the command to clean up          
        writetofile(remove_sh, sbatch(remove_sh,1,the_cwd) + [f'rm -rf {bamtmpdir} {splitsdir}\n', f'gzip ./{aligndir}/*.txt\n'], debug)
        ## Append the clean up command to file
        command_files.append((remove_sh,'afterok:',run_name,experi_mode,'clean','',remove_report))
    else: ## Otherwise do nothing
        pass
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      PIPELINE COMMAND SUBMISSION TO SLURM    
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      1) RESTARTING                          
    ## Call the command dataframe, remove previous logs if hard reset was passed 
    command_files, was_hard_reset = commandcontrol(command_files,hardreset,hic_pipeline,rerun,bwaix_jobid)
    ## 
    ##      2) SUBMITTING FASTP SPLITTING          
    ## Call the submission ftn for fastp commands 
    command_files, sub_sbatchs    = submitfastp(command_files,sub_sbatchs,partition,stamp,debugmode=debug)
    ## 
    ##      3) SUBMITTING BWA ALIGNMENTS           
    ## Call the submit bwa ftn
    command_files, sub_sbatchs    = submitbwa(command_files,sub_sbatchs,partition,stamp,bwa_runs,debugmode=debug)
    ## 
    ##      4) SUBMITTING PRE-FILTER COMMANDS
    ## Submit the pre commands and append the submitted jobs
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'pre','bwa',stamp,partition,bylast=True,debug=debug)
    ## 
    ##      5) SUBMITTING POST-FILTER COMMANDS      
    ## Submit the post commands and append the submitted jobs
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'post','pre',stamp,partition,bylast=True,debug=debug)
    ##
    ##      6) SUBMITTING Hi-C FILTER COMMANDS      
    ## Submit the hic filtering commands and append to the list of submitted jobs 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'filter','post',stamp,partition,bylast=True,debug=debug)
    ## 
    ##      7) SUBMITTING CONCAT COMMANDS         
    ## Submit the merge commands and append to sub 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'concat','filter',stamp,partition,group='Experiment' if postmerging else 'Sample',debug=debug)
    ## 
    ##      8) SUBMITTING SPLIT COMMANDS           
    ## Submit the chromosome splitting command
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'split','concat',stamp,partition,group='Experiment' if postmerging else 'Sample',debug=debug)
    ## 
    ##      9) SUBMITTING SORTTING COMMANDS         
    ## Submit the sort commands and append to sub 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'sort','split',stamp,partition,group='Experiment' if postmerging else 'Sample',debug=debug) 
    ## 
    ##      10) Hi-C file creation
    ## Call the juicer pre command for hic file creation if jarpath was passed 
    jar_batch = submitdependency(command_files,'juicerpre','sort',stamp,partition,group='Experiment' if postmerging else 'Sample',debug=debug) if jargiven else []
    ## Add to sub sbatchs 
    sub_sbatchs = sub_sbatchs + jar_batch
    ## Setnext step 
    next_step = 'juicerpre' if jargiven else 'sort'
    ## 
    ##      11) SUBMITTING COUNT COMMANDS 
    ## Submit the count command 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'count','sort',stamp,partition,group='Experiment',debug=debug) 
    ##
    ##      12) SUBMITTING TIME COMMANDS 
    ## Submit time stamp 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'timestamp',next_step,stamp,partition,group='Experiment',debug=debug) 
    ## 
    ##      13) CLEAN UP COMMANDS 
    ## Submit the clean up command if the flag was passed 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'clean','timestamp',stamp,partition,group='Experiment',debug=debug) 
    ## That is what we call a twelve step program :-)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    

    ##      SAVING OUT COMMAND FILE 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
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
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## End of file 