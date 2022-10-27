#!/usr/bin/env python

"""
Description: 
Removes all the gi numbers from fasta headers in old files downloaded from genbank
Usage: 
Run from folder enclosing the folder containing the suspect file(s).
Command line: fix_old_fasta_headers targetFolder
"""
#********************* MODULES**********************
from __future__ import division
from __future__ import print_function
import sys
import os
#******************* FUNCTION DEFINITIONS *******************



#********************** Main *************************

targetDir  = sys.argv[1]
os.chdir(targetDir)
Files = os.listdir('.')
for file in Files:
	if not file.startswith('.'):
		lineList = [] #list of all lines in the file
		INFILE = open(file, "rU")
		for line in INFILE:
			if line.startswith('>gi|'):
				temp = line.split(' ')
				temp2 = temp[0].split('|')
				AccNum = '>' + temp2[3]
				temp[0] = AccNum
				line = ' '.join(temp)
				lineList.append(line)
			else:
				lineList.append(line)
		INFILE.close()
		OUTFILE = open(file, 'w')
		for i in range(len(lineList)):
			OUTFILE.write(lineList[i])
		OUTFILE.close()
		
