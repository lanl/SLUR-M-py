#!/usr/bin/env python
## Import sys 
import sys 

## Bring in write to file 
from pysamtools import writetofile

## Set the output message and output file 
outmessage, outfile = ' '.join(sys.argv[1:-1]), sys.argv[-1]

## Write to file
writetofile(outfile,outmessage,debug=False)
## End of file 