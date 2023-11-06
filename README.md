# SLUR-M-py
SLUR(M)-PY; (pronounced slurpy) a slurm-pythonic pipeline for the alignment and analysis of paired-end sequencing data; from ATAC-seq to Hi-C
## Setting up the computing environment
Slurpy was developed using anaconda (python v 3.10.13). 
We recommend using conda to manage the python environment needed by slurpy.
Below are commands needed to set up the "bioenv" for running slurpy. 

```
## Make a python environment named bioenv 
conda create -n bioenv 

## Activate the environment
conda activate bioenv 

## Use conda to install needed mods
conda install numpy pandas matplotlib seaborn pysam samtools bwa macs2
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
### Envoking slurpy
The help menu (-h) of slurpy lists all the available arguments and default settings. Be sure to activate the conda environment "bioenv". 
```
## Activate our computing environment
conda activate bioenv 

## Call the help menu
./SLURPY/slurpy.py -h 

Calls the slurpy pipeline ( v 8.0.0 ) for alignment of pair-end reads to a reference genome. Example call: ./slurpy -e atac -r ../../REFERENCES/ENCODEREF/AUTOSOMES/GRCh38.fasta

options:
  -h, --help            show this help message and exit
  -r ./path/to/reference.bwaix, --refix ./path/to/reference.bwaix
                        Path to input reference bwa index used in analysis.
  -F 16, --fastp-splits 16
                        The number of splits to make for each pair of input fastq files (default: 16). Controls the total number of splits across the run.
  -B 16, --parallel-bwa 16
                        Number of parallel bwa alignments to run (default: 16). Controls the number of bwa jobs submitted at once to slurm.
  -P tb, --partition tb
                        The type of partition jobs formatted by slurpy run on (default: tb).
  -M chrM, --mtDNA chrM
                        Name of the mitochondrial contig (default: chrM).
  -Q 30, --map-threshold 30
                        Mapping quality threshold to filter alignments (default: 30).
  -f 16, --fastp-threads 16
                        The number of threads used in fastp to split input fastq files (default: 16). Note: must be an even multiple of the number of splits.
  -b 4, --bwa-threads 4
                        The number of threads used per bwa alignment on split input fastq files (default: 4).
  -s 4, --sam-threads 4
                        The number of threads used per bwa alignment on split input fastq files (default: 4).
  -n name, --run-name name
                        Run name used to name output files. Default behavior is to use the current parent directory.
  -c ./path/to/control.bam [./path/to/control.bam ...], --controls ./path/to/control.bam [./path/to/control.bam ...]
                        Path to control or input bam files used in ChIP-seq experiments.
  -g bp, --genome-size bp
                        Size of the genome being analyzed, used as parameter for macs2. Inputs can be integers in bp or two letter short hand, for e.g. hs for homo sapiens.
  -e type, --experiment type
                        The experiment type (default: wgs). Valid options include: atac, chip, wgs, hic.
  -R step, --rerun-from step
                        Step within the pipeline to re-run from. Options (in order) include: fastp, bwa, split, merge, mark, filter, macs2, count, clean.
  --restart             Flag to force the pipeline to reset and run from start.
  --runlocal            Disables sbatch submission and submits the script via bash to a local os.
  --debug               A flag to run in verbose mode, printing sbatch commands. Default behavior is false.
  --skip-dedup          Pass this flag to skip marking and removing duplicates. Default behavior is false (conduct duplicate marking).
  --broad               Flag to call broad peaks using the --broad-cutoff=0.1 setting in macs2. See macs2 callpeak --help for more details.
  --clean               If included will run clean up script at end of run. The default behavior is false, can be run after pipeline.

```
### For ATAC-seq experiments
To call the slurpy pipeline to analyze an ATAC-seq experiment run:

```
./SLURPY/slurpy.py -e atac -r /path/to/reference/file.fasta
```

### For ChIP-seq experiments 
To run slurpy to analyze a ChIP-seq experiment run:

```
./SLURPY/slurpy.py -e chip -r /path/to/reference/file.fasta -c /path/to/control/or/input.bam
```

## Dependencies
Slurpy utilizes SLURM and was developed under version 21.08.8-2. The suit of tools in samtools is also required with the minimum version of 1.15.1. 

### Additional linux core commands:
* cat 
* rm
* echo
