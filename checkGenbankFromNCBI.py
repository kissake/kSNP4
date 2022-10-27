#!/usr/bin/env python

"""
Description: searches headers.annotate_lit file for all accession numbers, then serches
genbank_from_NCBI.gbk files with each accession number.  Writes a list of those accession
numbers that are NOT in genbank_from_NCBI.gbk.
Run from within run folder created by kSNP3

Usage: check_genbank_form_NCBI
	returns a file named missing_accession_numbers.txt
"""
#********************* MODULES**********************
from __future__ import division
from __future__ import print_function
import sys
import os
#******************* FUNCTION DEFINITIONS *******************



#********************** Main *************************
AccList = [] # list of accession nmumbers
headersList = 'headers.annotate_list' # file of fasta headers
GBfile = 'genbank_from_NCBI.gbk'

#get alll the accession numbers
INFILE =open(headersList,'r')
for line in INFILE:
	temp = line.split(' ')
	accNum = temp[1]
	AccList.append(accNum)
INFILE.close()

#search genbank_from_NCBI.gbk
INFILE = open(GBfile, 'r')
data = INFILE.read().splitlines() #data is a list of lines ing GBfile
OUTFILE = open('missing_accession_numbers.txt', 'w')
for i in range(len(AccList)):
	accNum = AccList[i]
	present ='Missing'
	for i in range (len(data)):
		if accNum in data[i]:
			present = 'Present'
			print(accNum, '\t',present)
			break
	if present== 'Missing':
		print(accNum, '\t',present)
		string = accNum + '\n'
		OUTFILE.write(string)
INFILE.close()
OUTFILE.close()
