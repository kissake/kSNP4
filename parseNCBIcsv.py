#!/usr/bin/env python

"""
Description: Parses a .csv file from NCBI
Usage: parseNCBIcsv -i infileName -o outfilename -p prefix [-n N] where N is an integer 
		limiting the number of genomes listed.
		-L sets the level to I, incomplete genomes
Example1: parseNCBIcsv -i Ecoli.csv -o EcoliList.txt -p Eco -n 20
Example2: parseNCBIcsv -i Ecoli.csv -o EcoliAssemblies.txt -p Eco -n 20 -L
"""
#********************* MODULES**********************
import sys
import os
import re  #permits splitting on more than one charater at a time
import operator  #permits sorting multidimensional list
import random
import copy

#******************* FUNCTION DEFINITIONS *******************
def readInfile(infileName, level):
	#variables
	strain = '' #field 2
	replicons = '' #field 9
	assemblyID = '' #field 5
	numRep = 0 #field 11, number of replicons
	genomesList = []
	
	#open the infile
	#INFILE = open(infileName, 'r')
	#skip the first line 
	line = INFILE.readline()
	
	if level == 'C':
		for line in INFILE:
			line = line.rstrip()
			if line.find('Chromosome') ==-1 and line.find('Complete') ==-1:
				pass
			else:
				temp = line.split(',')
				genomesList.append([temp[2], temp[6], temp[9], int(temp[11])])
	elif level == 'I':
		for line in INFILE:
			line = line.rstrip()
			if line.find('Contig') ==-1 and line.find('Scaffold') ==-1:
				pass
			else:
				temp = line.split(',')
				genomesList.append([temp[2], temp[6], temp[5], int(temp[11])])
			
	return(genomesList)	
	
def parseGenomesList (genomesList):
	#variables
	numGenomes = len(genomesList)
	numRep = 0
	repStr = ''
	assemblyID = ''
	
	
	for i in range(numGenomes):
		#clean  up strain
		genomesList[i][0]=genomesList[i][0].replace(' ', '_')
		genomesList[i][0]=genomesList[i][0].replace('.', '_')
		genomesList[i][0]=genomesList[i][0].replace('\"', '')
		#clean up level
		genomesList[i][1]= genomesList[i][1].replace('\"', '')
		if level == 'C':
			#clean up replicons
			repStr = genomesList[i][2]
			numRep = genomesList[i][3]
			if numRep==1:
					temp = re.split("[:/]", repStr)
					genomesList[i][2] = temp[1]
					if genomesList[i][2].endswith('\"'):
						 genomesList[i][2] = genomesList[i][2].rstrip('\"')
			else:
				if repStr.find('/') != -1:
					pass
				else:
					repStr = repStr.replace(';', '/')
				
				temp = re.split("[:/]", repStr)
				k= len(temp)
				repStr = ''
				for j in range(k):
					if j%2==1:		
						repStr = repStr +' ' + (temp[j])
						if repStr.startswith(' '):
							repStr = repStr.lstrip()
						if repStr.endswith('\"'):
							repStr = repStr.rstrip('\"')			
				genomesList[i][2]= repStr		
		elif level == 'I':
			#clean up assemblyID
			assemblyID = genomesList[i][2]
			assemblyID = assemblyID.replace('\"', '')
			assemblyID = assemblyID.replace('A', 'F')
			genomesList[i][2] = assemblyID
	return genomesList
#********************** Main *************************
#variables
infileName = ''
outfileName = ''
level = 'C'  #C means only complete genomes, I means only incomplete genomes
numGenomesToGet = 0
prefix = ''
genomesList = [] #col0 = Strain, col1 = level c0l 2 = replicons, col 3 = numGenomes



#get command line arguments
for i in range(1, len(sys.argv)):
	if sys.argv[i] == '-i':
		infileName = sys.argv[i+1]
	if sys.argv[i] == '-o':
		outfileName = sys.argv[i+1]
	if sys.argv[i] == '-n':
		numGenomesToGet = int(sys.argv[i+1])
	if sys.argv[i] == '-p':
		prefix = sys.argv[i+1]
	if sys.argv[i] == '-L':
		level = 'I' 

#check to be sure numGenomesToGet is not >100
if numGenomesToGet > 100:
	sys.stdout.write('You may not set the -n option to >100\nQuitting parseNCBIcsv\nStart over with a lower value for -n\n.')
	sys.exit()
	
#open the infile
INFILE = open(infileName, 'r')

#read the infile
genomesList = readInfile(infileName, level)

#parse and clean up the genomes list
genomesList = parseGenomesList(genomesList)
numGenomes = len(genomesList)

# if -n has been specified randomize the order of genomes in genomesList and keep thee first N genomes
if numGenomesToGet > 0:
	random.shuffle(genomesList)
	newList = genomesList[0:numGenomesToGet]
	genomesList = copy.deepcopy(newList)
numGenomes = len(genomesList)

#sort genomesList by strain
genomesList.sort(key = operator.itemgetter(0))

#add prefix if one is selected
if prefix != '':
	for i in range (numGenomes):
		genomesList[i][0] = prefix + '_' + genomesList[i][0]



#write the outfile
OUTFILE = open(outfileName, 'w')
if level== 'C':
	OUTFILE.write('Strain\tAccNums\n')
elif level ==  'I':
	OUTFILE.write('Strain\tAssemblyID\n')
for i in range (numGenomes):
	OUTFILE.write('{0}\t{1}\n'.format(genomesList[i][0],genomesList[i][2]))



INFILE.close()
OUTFILE.close()
