#!/usr/bin/env python

"""
Description: 
Kchooser determines an optimum value for K then determines at that value of K
the fraction of kmers from the shortest sequence that are present in all of
the genomes.

The input file is the kSNP4 input file

Usage: Kchooser4 -in myfile.in [-u decimalFraction where the decimal fraction is the unique kmer limit]
Usage Example: Kchooser4 -in tinyData.in -u 0.985

By default chooses as the optimum value of K the odd number that is greater than the value
for which >0.99 of the kmers of the median-length sequence are unique.  The user can choose
a lower limit if desired.  See the kSNP4 documentation for an explanation.  To invoke the option enter
for instance: Kchooser4 -in myfile.fasta  -u 0.96.

Output is written to a file named Kchooser_infile.report, where infile is the name of the input file

This version uses the median length genomes at both steps
Example: Kchooser4 -in Ecoli96.in.  Output file is Kchooser4_Ecoli96.report

"""
#********************* MODULES**********************
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import time
import math
import subprocess
import operator

#******************* FUNCTION DEFINITIONS *******************
def readInfile(infileName):
	#variables
	GenomesList = []
	path = '' #path to a genome fasta file
	genomeLen = 0 #length of the genome in base pairs
	ID = '' #ID of the genome
	INFILE = open(infileName, "r")
	for line in INFILE:
		genomeSeq = ''
		line = line.rstrip()
		temp = line.split()
		path = temp[0]
		ID = temp[1]
		for seq_record in SeqIO.parse(path, 'fasta'):
			theSequence = seq_record.seq
			genomeSeq = genomeSeq + str(theSequence)
		genomeLen = len(genomeSeq)
		GenomesList.append([path, genomeLen, ID])
	return(GenomesList)

def findUniqueMerFraction(Kvalue, fractionUniqueKmers,theFile):
	#variables
	numUniqueKmers = 0
	kmers = []
	theNumber = str(10000000)
	
	subprocess.run(['jellyfish', 'count', '-C', '-o', 'output', '-s', theNumber, '-t', '3', theFile, '-m', str(Kvalue)]) 
	numMers = 0
	OUTFILE = open('jellyout.txt', 'w')
	subprocess.run(['jellyfish', 'dump', 'output_0', '-c'], stdout=OUTFILE) 
	OUTFILE.close()
	INFILE = open( 'jellyout.txt', 'r')
	for line in INFILE:
		line = line.rstrip()
		temp = line.split(' ')
		kmers.append([temp[0], temp[1]])
	numMers = len(kmers)
	for i in range(numMers):
		if kmers[i][1] == '1':
			numUniqueKmers += 1
	fractionUniqueKmers = numUniqueKmers/numMers
	INFILE.close()
	return (fractionUniqueKmers, numMers, kmers)

#*****************************************************	
#********************** Main *************************
#*****************************************************

#Variables
infileName  = ''
reportFileName = ''
fracUniqueKmerLimit = 0.99
GenomesList = [] #col0 is path to fasta genome file col1 is length of that genome
medianLength = 0 #median of genomes lengths
lenList = []
theSequence = '' #sequence of the median length genome
medIndex = 0 #index of the median length genome
thePath = '' #path to the median length genome
theID = '' #ID of median length genome
theGenomeLength = ''
Kvalue = 0 #the initial value of K, increments by 2 each round
KmerList = []
fractionUniqueKmers = 0
numGenomes = 0
Mers = 0
Fraction = 0
minGenomeID = ''
minGenomeLength = 0 #length of the shortest genome
minIndex = 0
numMers = 0
mers = []
sampleSize = 0
kmerList = [] #column 0 is a kmer, column 1 is the reverse complement of the kmer column 2 is the integer 1
revMer = ''
numTimes = 0 #the number of times the kmer is present in the median genome
coreKmers = 0
FCK = 0
previousFUK = 0 #the fractionUniqueKmers from the previous Kvalue
delta = 0 #difference between the current fractionUniqueKmers and previousFUK




startTime = time.time()

#get command line arguments
for i in range(1, len(sys.argv)):
	if sys.argv[i] == '-in':
		infileName = sys.argv[i+1]
	if sys.argv[i] == '-u':
		fracUniqueKmerLimit = float(sys.argv[i+1])

#read the input file
GenomesList = readInfile(infileName)
GenomesList.sort(key = operator.itemgetter(1)) #sorts GenomeList on column 1, the genome length
numGenomes = len(GenomesList)


#for i in range(numGenomes):
	#print(GenomesList[i][0], '\t', GenomesList[i][1])
	
#get the REPORTFILE name
temp = infileName.split('.')
reportFileName = 'Kchooser4_' + temp[0] +'.report'

#get the median length
print('Determining the median length of the genomes')
for i in range(len(GenomesList)):
	lenList.append(GenomesList[i][1])
medianLength = statistics.median(lenList)
medianLength = int(medianLength)

#get the length that is closest to the median length
for i in range (numGenomes):
	if lenList[i] - medianLength < 1: #
		closestLength  = lenList[i]
elapsedTime = time.time() - startTime
medianLength = closestLength


medIndex= lenList.index(medianLength)

thePath = GenomesList[medIndex][0]
theID = GenomesList[medIndex][2]
print('The median length genome is ',theID)
print('Its length is ', GenomesList[medIndex][1])

genomeSeq = ''
for seq_record in SeqIO.parse(thePath, 'fasta'):
			theSequence = seq_record.seq
			genomeSeq = genomeSeq + str(theSequence)
			

THEFILE = open('theFile', 'w')
THEFILE.write('>TheMedianFile\n')
THEFILE.write('{0}\n'.format(genomeSeq))
THEFILE.close()

elapsedTime = time.time() - startTime
print('Used {:.2f} seconds so far'.format(elapsedTime))

#get the initial value of K
theGenomeLength = GenomesList[medIndex][1]

Kvalue = (math.log(theGenomeLength) - math.log (0.1))/ math.log(4)
Kvalue = int(Kvalue)
if (Kvalue +1) %2 ==1:
	Kvalue = Kvalue+1
else:
	Kvalue = Kvalue+2
print('Initial value of K is ',Kvalue)

REPORTFILE = open(reportFileName, 'w')
REPORTFILE.write('Initial value of K is {0}\n'.format(Kvalue))

#find the optimum value of k
while fractionUniqueKmers <= fracUniqueKmerLimit:
	if previousFUK > 0:
		#print(previousFUK)
		(fractionUniqueKmers, numMers, mers) = findUniqueMerFraction(Kvalue, fractionUniqueKmers,'theFile')
		delta = fractionUniqueKmers-previousFUK	
		#print(' delta: ', delta)
		previousFUK = fractionUniqueKmers
		
	else:
		(fractionUniqueKmers, numMers, mers) = findUniqueMerFraction(Kvalue, fractionUniqueKmers,'theFile')
		previousFUK = fractionUniqueKmers
		delta = 1
		
		
	print('When k is ', Kvalue, ' ', fractionUniqueKmers, ' of the kmers from the median length sequence are unique.')
	REPORTFILE.write('When k is {0} {1} of the kmers from the median length sequence are unique.\n'.format(Kvalue,fractionUniqueKmers))
	if Kvalue < 31 and delta > 0.001:
		Kvalue = Kvalue + 2
	elif delta <= 0.001:
		break
	elif Kvalue > 31:
		print('WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!')
		print('K has exceeded the maximum value of 31 allowed by the kmer-counting program jellyfish')
		print('before the fraction of unique kmers reached the default value 0f 0.99.')
		print('You need to set a lower limit by adding the -u option to the command line.')
		print('Choose a lower limit, e.g. -u 0.985 for example. See the User Guide section on KChooser.')
		print('Terminating Kchooser now.')
		sys.exit()
if previousFUK> 0.99:
	Kvalue = Kvalue-2
print('The optimum value of k is ', Kvalue)
REPORTFILE.write('The optimum value of k is {0}\n\n'.format(Kvalue))


#write additional information
REPORTFILE.write('There were {0} genomes.\n'.format(numGenomes))
REPORTFILE.write('The median length genome was {0}\n'.format(theID))
REPORTFILE.write('Its length is {0}\n'.format(GenomesList[medIndex][1]))
elapsedTime = time.time() - startTime
print('\nUsed {0:.2f} seconds so far'.format(elapsedTime))
REPORTFILE.write('The time used to determine the optimum k was {0:.2f} seconds.\n'.format(elapsedTime))

#find FCK
print('Calculating the fraction of kmers that are core')
if numMers < 1000:
	sampleSize = numMers
else:
	sampleSize = 1000
print('sampleSize is ', sampleSize)

#put the first sampleSize kmers into kmerList
INFILE = open('jellyout.txt', 'r')
count = 0
for line in INFILE:
	if count < sampleSize:
		temp = line.split()
		mySeq = Seq(temp[0])
		forwardMer = str(mySeq)
		#print(temp[1])
		revMer = mySeq.reverse_complement()
		revMer = str(revMer)
		numTimes = int(temp[1])
		#print('\t', numTimes)
		if numTimes == 1:  #only count the unique kmers
			kmerList.append([forwardMer, revMer, 1])
			count += 1
		else:
			pass
	else:
		break
print('The number of kmers is ', len(kmerList))

INFILE.close()

#get the path and ID of the shortest genome
minGenomeID = GenomesList[0][2]
minGenomeLength = GenomesList[0][1]
print('The shortest genomes is ', minGenomeID, ' its length is ', GenomesList[0][1])
REPORTFILE.write('The shortest genomes is {0} its length is {1}\n'.format(minGenomeID, GenomesList[0][1]))

#Determine for the shortest genome whether each kmer in kmerList is in that genome sequence
thePath = GenomesList[0][0]
#genomeID = GenomesList[i][2]

theGenomeSequence = ''

for seq_record in SeqIO.parse(thePath, 'fasta'):
		theSequence = seq_record.seq
		theGenomeSequence = theGenomeSequence + str(theSequence)
		

for i in range(len(kmerList)):
	if theGenomeSequence.find(kmerList[i][0]) == -1 and theGenomeSequence.find(kmerList[i][1]) == -1:
		kmerList[i][2] = 0
	else:
		pass
elapsedTime = time.time() - startTime
print('Finished the shortest genome, {0}, at {1:.1f} seconds.'.format(minGenomeID, elapsedTime))

#calculate FCK
for i in range(len(kmerList)):
	coreKmers = coreKmers + kmerList[i][2]

print('There are {0} core kmers.'.format(coreKmers))
FCK = coreKmers/len(kmerList)
print('FCK = ', FCK)
REPORTFILE.write('\nFrom a sample of {0} unique kmers {1} are core kMers.\n'.format(len(kmerList), coreKmers))
REPORTFILE.write('FCK = {0:.3f}.\n'.format(FCK))

#clean up
os.remove('output_0')
os.remove('theFile')
os.remove('jellyout.txt')

elapsedTime = time.time() -startTime
REPORTFILE.write('The total time used was {0:.0f} seconds\n'.format(elapsedTime))
REPORTFILE.close()

