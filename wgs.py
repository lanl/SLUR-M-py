#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      wgs.py; Paired-End Sequencing Processing & Analysis Pipeline 
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



"""
## List slurpy command for VERO analysis
"""
./SLURPY/wgs.py -r ../Cryo_ATAC/GreenMVA_Ref/Chlorocebus_sabeus_mva.fasta -M NC_008066.1 

"""

##      MODULE LOADING and VARIABLE SETTING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in ftns and variables from defaults
from defaults import *
## Load in ftns from other libraries
from pysamtools import checksam, writetofile, txttobam, outnames
## Bring in ftns from atac-seq analysis 
from peaks import fastpeel, prepbwamem, splitcommand

## Set description   
basic_descr = "Processing pipeline for paired-end sequencing data."
## Set the hardset hic mode 
experi_mode = 'pe'
## Define the pipeline step help messages
R_help = "Step within the pipeline to re-run from. Options for basic, paired-end alignment include: %s"%(', '.join(basic_pipeline))
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
"""
This section is largly empty b/c we want line nubmers to (relatively) match up across our three main protocols. 
"""

## Jokes
"""
No functions were harmed in the making of this pipeline. - Sorry this feels like a good place to write out some jokes.


Q: What did one snowman say to the other?
A: Smells like carrots.

"""

## Thank the funding sources 
## 4DG funding 
"""
The work conducted within this pipeline was (first) funded by the Los Alamos National Laboratory Directed Research Grant 
(20210082DR) awarded to SRS and KYS.
"""

## BRAVE funding 
"""
This material is based upon work supported by the U.S. Department of Energy, Office of Science, 
through the Biological and Environmental Research (BER) and the Advanced Scientific Computing Research (ASCR) programs 
under contract number 89233218CNA000001 to Los Alamos National Laboratory (Triad National Security, LLC).
"""



































## All three python scripts, wgs.py (this script), atac.py, and hic.py share the same start line number below (106)
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is envoked by name 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = basic_descr)

    ## Add the required argument
    parser.add_argument("-r", "--refix",         dest="r", type=str,  required=True, help = r_help,  metavar = './path/to/reference.bwaix'       )
    
    ## Add the default arguments
    parser.add_argument("-F", "--fastp-splits",  dest="F", type=int,  required=False, help = F_help, metavar = splitsize,   default = splitsize  )
    parser.add_argument("-B", "--parallel-bwa",  dest="B", type=int,  required=False, help = B_help, metavar = parallelbwa, default = parallelbwa)
    parser.add_argument("-P", "--partition",     dest="P", type=str,  required=False, help = P_help, metavar = part,        default = part       ) 
    parser.add_argument("-M", "--mtDNA",         dest="M", type=str,  required=False, help = M_help, metavar = mito,        default = mito       )
    parser.add_argument("-R", "--rerun-from",    dest="R", type=str,  required=False, help = R_help, metavar = 'step',      default = None       )
    parser.add_argument("-q", "--fastq",         dest="q", type=str,  required=False, help = q_help, metavar = '.fastq.gz', default = fends      )

    ## Set number of threads 
    parser.add_argument("-f", "--fastp-threads", dest="f", type=int,  required=False, help = f_help, metavar = splitsize,   default = splitsize  )
    parser.add_argument("-b", "--bwa-threads",   dest="b", type=int,  required=False, help = b_help, metavar = bwathreads,  default = bwathreads )
    parser.add_argument("-t", "--sam-threads",   dest="t", type=int,  required=False, help = s_help, metavar = samthreads,  default = samthreads )

    ## Set values for run
    parser.add_argument("-n", "--run-name",      dest="n", type=str,  required=False, help = n_help, metavar = 'name',      default = None       )

    ## Set boolean flags 
    parser.add_argument("--restart",             dest="start",  help = restart_help,  action = 'store_true')
    parser.add_argument("--debug",               dest="debug",  help = debug_help,    action = 'store_true')
    parser.add_argument("--skip-fastp",          dest="sfast",  help = skipq_help,    action = 'store_true')
    parser.add_argument("--clean",               dest="clean",  help = clean_help,    action = 'store_true')
    parser.add_argument("--merge",               dest="merge",  help = merge_help,    action = 'store_true')

    



    
    """
    Adding space here to match the line numbers across scripts; if I didn't mention that already.
    """
    
    
    
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
    rerun           = inputs.R       ##     Setp to rerun pipeline from 
    fend            = inputs.q       ##     End of the input fastq files 
                                     ##     
    ## Set threads                   ##   
    fastp_threads   = inputs.f       ##     Number of fastp threads
    bwa_threads     = inputs.b       ##     Number of threads in bwa alignments
    sam_threads     = inputs.t       ##     Number of samtools threads 
                                     ##
    ## Set values for macs2          ##
    run_name        = inputs.n       ##     The name of the samples 
                                     ##   
    ## Set boolean vars              ##
    hardreset       = inputs.start   ##     Resetart the slurpy run, removing previous
    debug           = inputs.debug   ##     Run in debug mode 
    sfastp          = inputs.sfast   ##     Flag to skip fastp filtering 
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
    ## Initilizse list for sbatch
    sub_sbatchs = []











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
    ## Make the directories 
    makedirectories(grouped_dirs)
    ## Gather and remove the old command files, error, and output logs
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

    """
    Adding space here to hold line counts
    """





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
        bwaix_jobid = fakejobid if (debug and runlocal) else submitsbatch(intext)

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


    ##      SAMPLE NAME HANDELING
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
        ## Gather the fastp command and report 
        fastp_coms, fastp_report = fastpeel(r1,r2,'hic' if sfastp else experi_mode,fastp_threads,fastp_splits)
        ## Write the command to file
        writetofile(fastp_command_file, sbatch(None,fastp_threads,the_cwd) + fastp_coms, debug)
        ## Append command to file
        command_files.append((fastp_command_file,'singleton',sample_name,experi_mode,'fastp','',fastp_report))
        ## Format the split fastq files
        fastq_splits = [(f'{splitsdir}/{formatn(i+1)}.{basename(r1)}', f'{splitsdir}/{formatn(i+1)}.{basename(r2)}') for i in range(fastp_splits) ]
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      PARALLIZATION of ALIGNMENT 
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


            ##      BAM SPLITTING (post bwa) 
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## Format report name, commands, and file name 
            split_filename = f'{comsdir}/split.{basenobam(sbam)}.sh'
            ## Call the split ftn
            split_coms, split_report = splitcommand(sbam,mito,sam_threads)
            ## Write the split command
            writetofile(split_filename, sbatch(split_filename,sam_threads,the_cwd) + split_coms, debug)
            ## Apppend the seperating command
            command_files.append((split_filename,'afterok:',sample_name,experi_mode,'split','',split_report))
            ## Gather output bam name and txt names
            name_mapd_txt, name_plac_txt, name_mito_txt, name_unmp_txt, bedpe_txt = outnames(sbam,mito)    
            ## Remap the names to a bam file name 
            mapd_bam, plac_bam, mito_bam = txttobam(name_mapd_txt), txttobam(name_plac_txt), txttobam(name_mito_txt)
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


        ##      BAM MERGING 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Format commands to merge on mapped, placed, and mitochondrial 
        for bamend in ['mapped','placed',mito]:
            ## format the merged bam name and wild card
            tmp_inwild_bam = f'*.{bamend}.bam' if postmerging else f'*.{sample_name}.{bamend}.bam'
            ## Format the merge command 
            merge_coms, merge_report = mergebam(f'{sample_name}.{bamend}.bam', tmp_inwild_bam, sam_threads)
            ## Format the name of the merge command 
            merge_file = f'{comsdir}/concat.{sample_name}.{bamend}.sh'
            ## Write the merge command to file 
            writetofile(merge_file, sbatch(merge_file,sam_threads,the_cwd) + merge_coms, debug)
            ## Append the merge command
            command_files.append((merge_file,'afterok:',sample_name,experi_mode,'concat','',merge_report))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      MERGE UNMAPPED TXT 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Set the bam end file name
        bamend = 'unmapped'
        ## Format the infilename, we need this last bedpe for the hic
        newcatfile =  f'{aligndir}/{sample_name}.{bamend}.txt'
        ## Format the outfile 
        wildcard_filename = f'{bamtmpdir}/*.{bamend}.txt' if postmerging else f'{bamtmpdir}/*.{sample_name}*.{bamend}.txt'
        ## Format command to merge the unmapped txt file
        concat_coms, concat_report = pandacat(wildcard_filename,newcatfile)
        ## make concat file name
        concat_file = f'{comsdir}/concat.{sample_name}.{bamend}.sh'
        ## Write the concat command to file
        writetofile(concat_file, sbatch(concat_file,1,the_cwd) + concat_coms, debug)
        ## Append the concat command
        command_files.append((concat_file,'afterok:',sample_name,experi_mode,'concat','',concat_report))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## <- This indent marks all the following commands run post aligment and filtering of individual samples 
    

    ##      COUNTING BAM and TXT FILES      
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Format counting command, the count report, and the command
    counting_filename = f'{comsdir}/countbams.{run_name}.sh'                     ##      Name of the bash command file 
    count_report      = reportname(f'{aligndir}/bam.counts.bam','count')         ##      The name of the report 
    save_dist_name    = f'-S ./{diagdir}/{run_name}.bam.fragment.dist.png\n'     ##      Format the save file name
    ## Format the fragment histogram distribution for these samples 
    frag_calc_commands = f'{scriptsdir}/fragmentdist.py -b ./{aligndir}/*mapped*.bam ' + save_dist_name
    ## List the count commands 
    count_commands = [f'{scriptsdir}/countbams.py {run_name}\n', frag_calc_commands,  f'{scriptsdir}/myecho.py Finished counting bam files in {aligndir} dir. {count_report}\n']
    ## Wriet the coutn command to file
    writetofile(counting_filename,sbatch(counting_filename,1,the_cwd) + count_commands, debug)
    ## Append the counting command
    command_files.append((counting_filename,'afterok:',run_name,experi_mode,'count','',count_report))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SUBMITTING TIME-STAMP   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set the timesampe assocaited file names 
    timestamp_file   = f'{diagdir}/{run_name}.timestamp.{stamp}.txt'    ##     Name of the output file 
    timestampsh      = f'{comsdir}/time.stamp.sh'                       ##     Name of the .sh bash file 
    timestamp_report = reportname(run_name,f'timestamp.{stamp}')        ##     Name of the log to report to 
    ## Formath time stamp and echo commands 
    times_commands = [f'{scriptsdir}/endstamp.py {timestamp_file} {stamp}\n',  f'{scriptsdir}/myecho.py Finished SLURPY run of sample: {run_name}. {timestamp_report}\n']
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
    command_files, was_hard_reset = commandcontrol(command_files,hardreset,pipeline_steps,rerun,bwaix_jobid)
    ## 
    ##      2) SUBMITTING FASTP SPLITTING          
    ## Call the submission ftn for fastp commands 
    command_files, sub_sbatchs = submitfastp(command_files,sub_sbatchs,partition,stamp,debugmode=debug)
    ## 
    ##      3) SUBMITTING BWA ALIGNMENTS           
    ## Call the submit bwa ftn
    command_files, sub_sbatchs = submitbwa(command_files,sub_sbatchs,partition,stamp,bwa_runs,debugmode=debug)
    ## 
    ##      4) SUBMITTING SPLITTING COMMANDS
    ## Submit the bam split commands and append the submitted jobs
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'split','bwa',stamp,partition,bylast=True,debug=debug)
    ## 
    ##      5) SUMITTING THE CONCATONATION COMMANDS   
    ## Submit the concat commands and append the submitted jobs
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'concat','split',stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ##
    ##      6) SUBMITTING COUNTING COMMANDS      
    ## Submit the counting commands and append to the list of submitted jobs 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'count','concat',stamp,partition,group='Experiment',debug=debug)
    ##
    ##      7) SUBMITTING TIME COMMANDS 
    ## Submit time stamp 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'timestamp','count',stamp,partition,group='Experiment',debug=debug)
    ## 
    ##      8) CLEAN UP COMMANDS 
    ## Submit the clean up command if the flag was passed 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'clean','timestamp',stamp,partition,group='Experiment',debug=debug)
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