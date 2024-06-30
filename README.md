# SLURPY
[SLUR(M)-py: A SLURM Powered Pythonic Pipeline for Parallel Processing of 3D (Epi)genomic Profiles](https://www.biorxiv.org/content/10.1101/2024.05.18.594827v1)
## Setting up the computing environment
Slurpy was developed using anaconda (python v 3.10.13). 
We recommend using conda to manage the python environment needed by slurpy.
Below are commands needed to set up the "bioenv" for running slurpy. 

```
## Make a python environment named bioenv 
conda create -n bioenv 

## Activate the environment
conda activate bioenv 
```

After making a new conda environment install needed packages.

```
## Use conda to install needed mods
conda install numpy pandas matplotlib seaborn pysam samtools biopython dask samblaster bwa macs2
```

If the above installation command hangs, we recommend removing macs2 from the list of libraries and trying again. Then installing macs2 via pip.

```
pip install macs2
```
A full list of the python libraries and their versions used to develop slurpy are listed within [python.dependencies.txt](https://github.com/SLUR-m-Py/SLURPY/blob/main/python.dependencies.txt).

## Installation
Downloading the repository as a .zip archive is easiest. For developers a simple clone command with git works too:

```
git clone https://github.com/SLUR-m-Py/SLURPY.git
```

Once slurpy is downloaded (and expanded), change the current directory to the local SLURPY directory and modify the python files as executables. 

```
cd ./SLURPY
chmod +x *.py 
```

### Checking the python environment 
Once within the slrupy directory, run the environment checking script, [modcheck.py](https://github.com/SLUR-m-Py/SLURPY/blob/main/modcheck.py).

```
## Activate the computing environment
conda activate bioenv 
python modcheck.py
```

If the environment was created successfully the script will run to completion and print the following:

```
INFO: Modules loaded correctly.
INFO: We are ready to run slurpy!
```

## Running slurpy
### Setting up a project directory
Currently the ATAC- and ChIP-seq protocols are fully functional. To run slurpy, change the current directory to the target project (in the example below, the project is named 2501_001) and soft link to the path of the slurpy executables.

```
cd /path/to/project/directory/2501_001
ln -s /path/to/SLURPY
```

Within the project directory ensure the paired fastqs are listed or linked within a subdirectory 
“./fastqs”. 

```
ls -l ./fastqs/*.fastq.gz
```
### Envoking slurpy for Hi-C processing
The help menu (-h) of protocols within slurpy lists all the available arguments and default settings. Be sure to activate the conda environment "bioenv". 
```
## Activate our computing environment
conda activate bioenv 

## Call the help menu for hic.py 
./SLURPY/hic.py -h
usage: hic.py [-h] -r ./path/to/reference.bwaix [-F 64] [-B 64] [-P tb] [-M chrM] [-Q 30] [-R step] [-q .fastq.gz] [-f 8] [-b 4] [-t 4] [-n name] [-E bp] [-C bp] [-L MboI] [-D n] [-Z n] [-G ./path/to/list.tsv] [-J ./path/to/juicer.jar]
              [-S 25000, 10000, ... [25000, 10000, ... ...]] [--restart] [--debug] [--skip-dedup] [--clean] [--merge]

Processing and analysis pipeline for paired-end sequencing data from Hi-C experiments.

optional arguments:
  -h, --help            show this help message and exit
  -r ./path/to/reference.bwaix, --refix ./path/to/reference.bwaix
                        Path to input reference bwa index used in analysis.
  -F 64, --fastp-splits 64
                        The number of splits to make for each pair of input fastq files (default: 64). Controls the total number of splits across the run.
  -B 64, --parallel-bwa 64
                        Number of parallel bwa alignments to run (default: 64). Controls the number of bwa jobs submitted at once to slurm.
  -P tb, --partition tb
                        The type of partition jobs formatted by slurpy run on (default: tb).
  -M chrM, --mtDNA chrM
                        Name of the mitochondrial contig (default: chrM).
  -Q 30, --map-threshold 30
                        Mapping quality threshold to filter alignments (default: 30).
  -R step, --rerun-from step
                        Step within the pipeline to re-run from. Options for Hi-C analysis include: fastp, bwa, pre, post, filter, concat, split, sort, count, clean
  -q .fastq.gz, --fastq .fastq.gz
                        The file extension of input fastq files (default: .fastq.gz)
  -f 8, --fastp-threads 8
                        The number of threads used in fastp to split input fastq files (default: 8). Note: must be an even multiple of the number of splits.
  -b 4, --bwa-threads 4
                        The number of threads used per bwa alignment on split input fastq files (default: 4).
  -t 4, --dask-threads 4
                        The number of threads used in calls to functions and calculations with pandas and dask dataframe(s) (default: 4).
  -n name, --run-name name
                        Run name used to name output files. Default behavior is to use the current parent directory.
  -E bp, --error-distance bp
                        Linear genomic distance to parse left and right oriented, intra-chromosomal Hi-C pairs for missing restriciton site(s). Passing zero (0) will skip this check (default: 10000 bp).
  -C bp, --self-circle bp
                        Linear genomic distance to check outward facing, intra-chromosomal Hi-C contacts for self-circle artifacts. Passing zero (0) will skip this check (default: 30000 bp).
  -L MboI, --library MboI
                        The name of the restriction site enzyme (or library prep) used in Hi-C sample creation. Options include Arima, MboI, DpnII, Sau3AI, and HindIII (default: Arima). Passing none (i.e. Dovetail) is also allowed, but checks for
                        restriction sites and dangling ends will be skipped.
  -D n, --mindist n     A filter on the minimum allowed distance (in bp) between reads (within a pair) that make up an intra-chromosomal Hi-C contact. Default behaviour is none (i.e. default: 0).
  -Z n, --chunksize n   Number of rows (default: 50000) loaded into pandas at a time. WARNING: while increasing could speed up pipeline it could also cause memeory issues.
  -G ./path/to/list.tsv, --genomelist ./path/to/list.tsv
                        Path to list of chromosomes (by name) to include in final Hi-C analysis. Must be a tab seperated tsv or bed, comma seperated csv, or space seperated txt file with no header.
  -J ./path/to/juicer.jar, --jar-path ./path/to/juicer.jar
                        Path to juicer jar file for juicer pre command. Required for .hic file creation.
  -S 25000, 10000, ... [25000, 10000, ... ...], --bin-sizes 25000, 10000, ... [25000, 10000, ... ...]
                        Chromosome resolution (i.e. bin sizes) for .hic files. Default: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000
  --restart             Flag to force the pipeline to reset and run from start.
  --debug               A flag to run in verbose mode, printing sbatch commands. Default behavior is false.
  --skip-dedup          Pass this flag to skip marking and removing duplicates. Default behavior is false (conduct duplicate marking).
  --clean               If included will run clean up script at end of run. The default behavior is false, can be run after pipeline.
  --merge               Passing this flag will merge across all pairs of fastqs for final output.
```
### For ATAC-seq experiments
To call the peaks.py script within the slurpy pipeline to analyze an ATAC-seq experiment run:

```
./SLURPY/peaks.py -r /path/to/reference/file.fasta
```

### For ChIP-seq experiments 
To run slurpy to analyze a ChIP-seq experiment run:

```
./SLURPY/peaks.py -r /path/to/reference/file.fasta -c /path/to/control/or/input.bam
```

## Dependencies
Slurpy utilizes SLURM and was developed under version 21.08.8-2. The suit of tools in samtools is also required with the minimum version of 1.15.1. 

### Additional linux core commands:
* cat 
* rm
* echo
