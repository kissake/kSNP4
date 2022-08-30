#!/usr/bin/env python

"""
Description: Downloads a set of fasta genome sequences. Infile is an output file from parseNCBIcsv.py

Usage: getFastaGenomes -i infileName -o outDir -e email
outDir is the output directory to hold all downloaded genomes. It must be in same 
directory from which script is run.
email is your email address and is required
Example: getFastaGenomes -i ShortEcoliList.txt -o GenomeFiles -e myName@gmail.com
"""
#********************* MODULES**********************
import sys
import os
from Bio import Entrez
import time
#******************* FUNCTION DEFINITIONS *******************
def readInfile(infileName):
	#variables
	genomesList = []
	ID = ''
	accNums = ''
	
	#open the infile
	INFILE = open(infileName, 'r')
	
	#skip the first line 
	line = INFILE.readline()
	for line in INFILE:
		line = line.rstrip()
		temp = line.split('\t')
		ID = temp[0]
		accNums = temp[1]
		genomesList.append([ID,accNums])		
	INFILE.close()
	return genomesList

def getGenome(ID, accList, outDir):
	#variables
	handle='' # handle for Entrez.efetch
	numRep = 0 #number of replicons in the genome
	fileName = ID +'.fna'
	OUTFILE = open(fileName, 'w')
	#print(ID)
	#print(accList)
	#split accList into its components
	temp = accList.split(' ')
	numRep = len(temp)
	startTime = time.time()
	for i in range(numRep):
		handle = Entrez.efetch(db="nuccore", ID = temp[i], rettype = "fasta", retmode = "text")
		OUTFILE.write("{0}\n\n".format(handle.read()))
	elapsedTime = time.time() - startTime
	sys.stdout.write("Time to download {0} was {1:.2f} seconds.\n".format(fileName, elapsedTime))
	OUTFILE.close()
#********************** Main *************************
#variables
infileName = '' #name of the input file
outDir = ''  #name of the output directory
email = ''
genomesList = [] #col0= genome name, col1 = accession numbers separated by space
numGenomes = 0
home = '' #the directory from which the program is reun


#get command line arguments
for i in range(1, len(sys.argv)):
	if sys.argv[i] == '-i':
		infileName = sys.argv[i+1]
	if sys.argv[i] == '-o':
		outDir = sys.argv[i+1]
	if sys.argv[i] == '-e':
		email = sys.argv[i+1]

#check that all required arguments have been entered on the command line
if infileName == '':
	sys.stdout.write('You must enter the infile name on the command line as -i infileName.\nQuitting the program now.\n')
	sys.exit()
if outDir == '':
	sys.stdout.write('You must enter the output directory name on the command line as -o outputDirectoryName.\nQuitting the program now.\n')
	sys.exit()
if email == '':
	sys.stdout.write('You must enter your email address on the command line as -e myself@someplace.com.\nQuitting the program now.\n')
	sys.exit()
	
	
home = os.getcwd()

#read the infile
genomesList = readInfile(infileName)
numGenomes = len(genomesList)

#send the email address to NCBI
Entrez.email = email

#move to the outdir
os.chdir(outDir)

#download the fasta records from Entrez
for i in range(numGenomes): #numGenomes
	getGenome(genomesList[i][0], genomesList[i][1], outDir)
