#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## -------------------------------- v 0.0.0 --------------------------------- ## 
##      HARDRESET: Completely deletes directories and files made by SLURPY    ##
## -------------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## -------------------------------------------------------------------- ##
##      MODULE LOADING 
## Load in shutil
import shutil       
## Load in direcotries
from defaults import diagdir, aligndir, bamtmpdir, debugdir, splitsdir, comsdir
## -------------------------------------------------------------------- ##

## -------------------------------------------------------------------- ##
##      RESETING 
## Remove each dir
[shutil.rmtree(d,ignore_errors=True) for d in [diagdir, aligndir, bamtmpdir, debugdir, splitsdir, comsdir]]
## -------------------------------------------------------------------- ##
## End fo file 