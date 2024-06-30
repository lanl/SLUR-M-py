#!/usr/bin/env python
## List of needed python libraries
"""
numpy
pandas
matplotlib
seaborn 
"""
## List of standard libraries with python 3
"""
glob
os
subprocess
argparse
time
sys
json
multiprocessing
datetime
biopython
dask
"""
## Needed stand alones 
"""
samtools
bwa
macs2
"""
## Conda install command 
"""
conda install numpy pandas matplotlib seaborn ipython pysam samtools macs2 bwa 
"""
## If the conda install dosn't work for macs2, try:
"""
pip install macs2
"""
## Load in date and time
from datetime import datetime

## Load in glob
from glob import glob 

## Load in SeqIO
from Bio import SeqIO

## Bring in basename, get file size, and path-exists 
from os.path import isfile, basename, exists, getsize 

## Bring in make dirs
from os import makedirs, remove, getcwd

## bring in matplot lib 
from matplotlib import pyplot as plt

## Load pandas, numpy and seaborn 
import pandas as pd, numpy as np, seaborn as sns 

## Load in system, json, time, argparese and subprocess 
import sys, json, time, argparse, subprocess

## Loadin checksam
from pysamtools import checksam

## Bring in dask 
import dask.dataframe as dd 

## Check samtools
assert checksam(), 'ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again.\n'

## Print to screen that this ran
print("INFO: Modules loaded correctly.\nINFO: We are ready to run slurpy!")
## End of file 