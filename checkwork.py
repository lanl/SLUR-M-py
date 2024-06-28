#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## ----------------------------- v 0.0.0 ------------------------------ ## 
##          CHECKWORK: Checks error logs from runs of SLURPY            ##
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
## -------------------------------------------------------------------- ##
##      MODULE LOADING 
## Bring in sorted glob 
from defaults import sortglob, debugdir
## Bring in get size from os 
from os.path import getsize  
## bring in numpy 
import numpy as np 
## -------------------------------------------------------------------- ##

## -------------------------------------------------------------------- ##
##      ANALYSIS of ERROR LOGS
## Bring in the error logs 
all_error_logs = sortglob(f'./{debugdir}/*.err')
## Filter the error logs for non bwa and fastp logs
error_logs =  np.array([k for k in all_error_logs if ('bwa' not in k) and ('fastp' not in k)])
## Gather the sizes 
error_sizes = np.array([getsize(k) for k in error_logs])
## -------------------------------------------------------------------- ##

## -------------------------------------------------------------------- ##
##      CHECK ERROR STATUS
## Check if we have non-zeros
if np.sum(error_sizes):
    ## Gather the logs with erros 
    has_errors = error_logs[(error_sizes > 0)]

    ## Print the number of errors
    print(f'WARNING: Errors were detected in {len(has_errors)} log(s)!')
    print("WARNING: Check the following log(s):\n")
    ## Print the logs with erros
    [print("\t%s"%f) for f in has_errors]

else: ## Otherwise
    print("INFO: Good news everyone! No errors were detected in this run.\n\t:-)\n\t<3")
## -------------------------------------------------------------------- ##
## End of file 