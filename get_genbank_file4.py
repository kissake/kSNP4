#!/usr/bin/env python

"""
Description: 
AccNum.txt is a file of accession numbers

Usage: 
GenbankDownload pathToThegenomeFasta file
"""
#********************* MODULES**********************
from __future__ import division
from __future__ import print_function
import sys
import os
import os.path
import sys
import re #regular expressions
import time
from Bio import Entrez
#******************* FUNCTION DEFINITIONS *******************
def extractAccNum(line):
	temp = line.split()
	accnum = temp[0]
	accnum = accnum.lstrip('>')	
	AccNumList.append(accnum)
	print(accnum)
	return(accnum)

#********************** Main *************************
#variables
fileExist = ''
pathToFastaGenome = sys.argv[1]
accnum = ''
GenomeFile = ''
seqData = False  # et to True iuf line is sequence data
AccNumList = []

startTime = time.time()

#does the file genbank_from_NCBI.gb esist?
fileExist = os.path.isfile('./genbank_from_NCBI.gbk') #'genbank_from_NCBI.gbk'
#print(fileExist)
if fileExist == True: # if the file does  exist do nothing
	pass
elif fileExist == False: #if the file does not exist create it
	OUTFILE = open('genbank_from_NCBI.gbk', 'w')
	#OUTFILE.write('\n')
	OUTFILE.close()
	
#print(pathToFastaGenome)

#extract the accession number from the sequence file header
INFILE = open(pathToFastaGenome, 'rU')

for line in INFILE:
	if line.startswith('>'):
		accnum = extractAccNum(line)
INFILE.close


# download the accnum file
#handle = Entrez.efetch(db="nuccore", id=accnum, retmax="100", rettype="gb", retmode="text") #, retmode="text"
#text = handle.read() #Entrez.read(handle)

#temp = text.split('\n')

OUTFILE = open('genbank_from_NCBI.gbk', 'a',)

for i in range(len(AccNumList)): #for each accession number
	accnum = AccNumList[i]
	#find the NCBI ID for the accnum
	Entrez.email = 'barryghall@gmail.com'
	handle = Entrez.esearch(db="nuccore", term=accnum, retmax=100)
	records = Entrez.read(handle) ##records is a dict
	# print(records, '\n\n')
	#identifiers = records['IdList']
	ID = records['IdList']
	if len(ID) > 0:
		ID = str(ID[0])
		#print("The ID is ", ID)
		try:
			handle = Entrez.efetch(db="nuccore", id=accnum, retmax="100", rettype="gb", retmode="text") #, retmode="text"
		except:
			continue
		text = handle.read() #Entrez.read(handle). text is a long string with the contents of downloaded file
		temp = text.split('\n')	
		for j in range(len(temp)):	
			if temp[j].startswith("ORIGIN"):
				OUTFILE.write('{0}\n'.format(temp[j]))
				OUTFILE.write('//\n\n')		
				break
			else:
				 # download the accnum file
				 OUTFILE.write('{0}\n'.format(temp[j]))

		OUTFILE.write('//\n\n')
	time.sleep(1)

OUTFILE.close()
endTime = time.time()
elapsedTime = endTime-startTime
print("Used ",elapsedTime, "seconds.")

