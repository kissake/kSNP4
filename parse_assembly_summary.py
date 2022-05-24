#!/usr/bin/env python3

"""Description: 

Usage: parse_assembly_summary infileName outfileName

"""
#********************* MODULES**********************
from __future__ import division
from __future__ import print_function
import sys
import os
#******************* FUNCTION DEFINITIONS *******************



#********************** Main *************************
infile = sys.argv[1] #infile name
outfile = sys.argv[2]

#path is element 19, genome ID is element 7, type is element 11
#read the infile
INFILE = open(infile, "rU")
OUTFILE = open(outfile, "w")

for line in INFILE:
	if line.startswith('#'):
		pass
	else:
		line = line.rstrip()
		temp = line.split('\t')
		ID = temp[7]
		path = temp[19]
		type = temp[11]
		OUTFILE.write("{0}\t{1}\t{2}\n".format(ID, type, path))
INFILE.close()

OUTFILE.close()

