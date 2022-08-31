#!/usr/bin/env python3

"""
Description: Parses the the genbank_from_NCBI.gbk and SNPs_all files written by kSNP3
to generate an annotations file

Usage: python ParAnn.py 
"""
#********************* MODULES**********************
from __future__ import division
from __future__ import print_function
import sys
import os
import re
import time
import resource #for determining peak memory used
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna  ## Removed from Bio module per: https://biopython.org/wiki/Alphabet
#******************* FUNCTION DEFINITIONS *******************
def findAA(locusSeq, SNPallele, SNPpos, start, end, strand, complement):
	temp = locusSeq.split('.')
	flankLen = len(temp[0])
	theSeq = temp[0] + SNPallele + temp[1]
	
	#find the reading frame	
	if strand == 'R':
		readingFrame = ((end+1) - SNPpos)%3
	elif strand == 'F':
		readingFrame = (SNPpos - (start-1))%3  #this calculation is correct
	
	if readingFrame == 0:
		codon = theSeq[flankLen-2:flankLen+1]
		peptideSeq = theSeq[flankLen-5:flankLen+7]
	elif readingFrame == 1:
		codon = theSeq[flankLen:flankLen+3]
		peptideSeq = theSeq[flankLen-3:flankLen+9]
	elif readingFrame == 2:
		codon = theSeq[flankLen-1:flankLen+2]
		peptideSeq = theSeq[flankLen-4:flankLen+8]
	codon = Seq(codon)
	peptideSeq = Seq(peptideSeq)
	
	RC = 'F' #don't reverse complement the codon or the peptideSeq
	if (strand == 'R' and complement == 'F' ):  # or (strand == 'F' and complement == 'T')
		codon = codon.reverse_complement()
		peptideSeq = peptideSeq.reverse_complement()
		RC='T'
	elif strand == 'F' and complement == 'T':
		codon = codon.reverse_complement()
		peptideSeq = peptideSeq.reverse_complement()
		RC='T'
	peptide = peptideSeq.translate()	
	aminoAcid = codon.translate()
	return(str(peptide),str(aminoAcid), str(codon), readingFrame, RC)

def getNewCodon(base,codon, readingframe, RC):
	if RC == 'T':
		base = Seq(base) #base is now a seq object
		base = base.reverse_complement()
		base = str(base)
	codonList = list(codon)
	if readingFrame == 0 and RC == 'F':
		codonList[2] = base
		newCodon = ''.join(codonList)
	elif readingFrame == 1 and RC == 'F':
		codonList[0] = base
		newCodon = ''.join(codonList)
	elif readingFrame == 2 and RC == 'F':
		codonList[1] = base
		newCodon = ''.join(codonList)
	elif readingFrame == 0 and RC == 'T':
		codonList[0] = base
		newCodon = ''.join(codonList)
	elif readingFrame == 1 and RC == 'T':
		codonList[2] = base
		newCodon = ''.join(codonList)
	elif readingFrame == 2 and RC == 'T':
		codonList[1] = base
		newCodon = ''.join(codonList)
	newCodon = Seq(newCodon)
	aa = newCodon.translate()

	return(str(newCodon), str(aa))
#********************** Main *************************
outfileName = 'SNPs_all_annotated'
summaryFile = 'Annotation_summary'
genbankFile = 'genbank_from_NCBI.gbk'
SNPs_allFile = 'SNPs_all'
RefGenomesDict = {} #key is accession number value is annoymous array of gene information
				#the annonymous array is [col0 = geneType, col1 = complement, col2 = start, col3 = end, col4 = product]
startTime = time.time()
endTime = 0
usedTime = 0
AccNumm = '' #the accession number of the genbank file currently being parsed
geneType = '' #either CDS or RNA
complement = '' # either T or F
start = 0 #start position of a gene
end = 0 #end position of a gene
product = '' #the product of a gene
SNPnum = 0  #number of the locus
locusSeq = '' #sequence of the SNP locus with central base indicated by dot
SNPallele = '' #the central base; i.e the SNP
GenomeID = "" # ID for the genome used by kSNP3
SNPpos = ''  #position of the SNP in the genome if in a reference genome,
strand = '' #strand on which the SNP locus occurs, F or R
SNPinfo = [] # two dimensional list col0=SNPnum, col1 = locusSeq, col2= SNPallele, col3 = GenomeID, col4=SNPpos, col 5 = strand, col6 = AccNum
infoList = []



#Get accession numbers and initialize
INFILE = open(genbankFile, "r")
for line in INFILE:
	line = line.rstrip()
	if line.startswith('VERSION'):
		temp = re.split(r" +", line)
		RefGenomesDict[temp[1]] = []
INFILE.close()		

INFILE = open(genbankFile, "r")
for line in INFILE:
	line = line.rstrip()
	if line.startswith('VERSION'):
		temp = re.split(r" +", line)
		AccNum = temp[1]
	elif re.search (r"CDS  ", line) or re.search (r"RNA  ", line) and re.search(r"\.\.", line): #if a coding sequence		
		temp = re.split(r" +", line)
		geneType = temp[1]
		Range = temp[2]
		if not re.search(r"[<>]",Range) and  not re.search(r"join",Range):
			if re.search(r"complement", Range):
				complement = 'T'
				temp =re.split(r"[()]", Range)
				#print(temp[1])
				temp2 = re.split(r"\.\.",temp[1])
				start = int(temp2[0])
				end = int(temp2[1])
			else:
				complement = 'F'
				temp = re.split(r"\.\.", Range)
				start = int(temp[0])
				end = int(temp[1])		
			while not re.search(r"/product", line):
				line = INFILE.readline()
				line = line.rstrip()
			temp = line.split("\"")
			product = temp[1]
			RefGenomesDict[AccNum].append([geneType, complement, start, end, product])
INFILE.close()


endTime = time.time()
usedTime = endTime - startTime

######################### parse the SNPs_all file and put info into SNPinfo list##########
INFILE = open(SNPs_allFile, "r")
line = INFILE.readline()  #skips the first line which is a blank
for line in INFILE:
	line = line.rstrip()
	if len(line) >0:  #if this is not a blank line
		temp = line.split('\t')
		SNPnum = int(temp[0])
		locusSeq = temp[1]
		SNPallele = temp[2]
		SNPpos = temp[3]
		GenomeID = temp[4]
		if len(temp) >5:
			AccNum = temp[5]
		else:
			AccNum = 'NoAccNum'
		if SNPpos == 'x':
			SNPinfo.append([SNPnum, locusSeq, SNPallele, GenomeID, SNPpos, '---',AccNum])
		else:
			temp2 = SNPpos.split(' ')
			SNPpos = int(round(float(temp2[0])))
			strand = temp2[1]
			SNPinfo.append([SNPnum, locusSeq, SNPallele, GenomeID, SNPpos, strand, AccNum])
INFILE.close()		
lastSNP = SNPinfo[-1][0]
		
#########################################
## now handle all that information and write the output file
OUTFILE = open(outfileName, "w")
SUMMARY = open(summaryFile, "w")
currentSNP = 0
codon = ''
readingFrame = -1
for info in SNPinfo:
	#if currentSNP >2590:
		#sys.exit()
	if 'infoList' in locals() and len(infoList) >0 and info[0] != currentSNP:

#****** start block 1 ******
		for i in range(len(infoList)):
			if InAnnotatedGenome == 'F':
				OUTFILE.write('{0}\tSNPallele: {1}\tGenomeID: {2}\tNot in annotated genome\n'.format(infoList[i][0], infoList[i][1], infoList[i][3]))
			else:
				if InAnnotatedRegion == 'F':
					OUTFILE.write('{0}\tSNPallele: {1}\tGenomeID: {2}\tNot in annotated region\n'.format(infoList[i][0], infoList[i][1], infoList[i][3]))
				else:
					if infoList[i][5]== 'rRNA' or infoList[i][5]== 'tRNA':
						OUTFILE.write('{0}\tSNP allele: {1}\tSNP position: {2}\tGenome ID: {3}, AccNum: {4}\t {5}\t{6}\n'.format(infoList[i][0],infoList[i][1],infoList[i][2], infoList[i][3], infoList[i][4], geneType,infoList[i][9]))
					else:	#if in a coding sequence
						OUTFILE.write('{0}\tSNP allele: {1}\t'.format(infoList[i][0], infoList[i][1])) #SNPnum and SNPallele
						OUTFILE.write('SNP position: {0}\t'.format(infoList[i][2])) #SNPpos						
						OUTFILE.write('Genome ID: {0}\tAccNum: {1}\t'.format(infoList[i][3],infoList[i][4] ))  #genomeID & accession number
						if infoList[i][6] != '?':
							OUTFILE.write('codon: {0}\t'.format(infoList[i][6]))  #codon
						else:
							if len(codonInfo) >1:
								thisCodon = codonInfo[infoList[i][1]]
								infoList[i][6] = thisCodon
								OUTFILE.write('codon: {0}\t'.format(infoList[i][6]))
							else:
								(codon, aa) = getNewCodon(infoList[i][1], codon, readingFrame, RC)
								infoList[i][6] = codon
								infoList[i][7] = aa
								OUTFILE.write('codon: {0}\t'.format(infoList[i][6]))
						
						if infoList[i][7] != '?':
							OUTFILE.write('amno acid: {0}\t'.format(infoList[i][7]))  #amino acid
						else:
							OUTFILE.write('amno acid: {0}\t'.format(infoList[i][7]))  #amino acid
						OUTFILE.write('peptide: {0}\tgene product: {1}\n'.format(infoList[i][8], product))  #codon
		#now report SNPnum, alternative SNP alleles, alternative codons, alternative amino acids, synonymnous or non synonymous and product				
		SNPid = infoList[0][0]
		#get the alternative SNPs, amino acids and codons
		altSNPs = []
		altCodons = []
		altAA = []
		for i in range(len(infoList)):
			if infoList[i][1] not in altSNPs:
				altSNPs.append(infoList[i][1])
			if infoList[i][6] != '?':
				if infoList[i][6] not in altCodons:
					altCodons.append(infoList[i][6])
			if infoList[i][7] != '?':	
				if infoList[i][7] not in altAA:
					altAA.append(infoList[i][7])
		#make the alternates into strings
		SNPstring = ''
		for i in range(len(altSNPs)):
			SNPstring =SNPstring + altSNPs[i] + '_'
		SNPstring = SNPstring[:-1]
		
		
		if len(altCodons)>0:
			codonString = ''
			for i in range(len(altCodons)):
				codonString = codonString + altCodons[i] + '_'
			codonString = codonString[:-1]
		
		if len(altAA)>0:
			aaString = ''
			for i in range(len(altAA)):
				aaString = aaString + altAA[i] + '_'
			aaString = aaString[:-1]
		
		SUMMARY.write('{0}\tSNPs: {1}\t'.format(SNPid,SNPstring))
		if InAnnotatedGenome == 'F':
			SUMMARY.write('Not in annotated genome\n')
		elif InAnnotatedGenome == 'T' and InAnnotatedRegion == 'F':
			SUMMARY.write('Not in annotated region of annotated genome.\n')
		elif InAnnotatedGenome == 'T' and InAnnotatedRegion == 'T':
			if geneType =='rRNA' or geneType == 'tRNA':
				SUMMARY.write('RNA\tGene product: {0}.\n'.format(product))
			else:
				if len(aaString) >1:
					SUMMARY.write('Codons: {0}\tAmino acids: {1}\tNonSynonymous\tGene product: {2}.\n'.format(codonString, aaString, product))
				else:
					SUMMARY.write('Codons: {0}\tAmino acid: {1}\tSynonymous\tGene product: {2}.\n'.format(codonString, aaString, product))
		##############
		OUTFILE.write('\n')
		currentSNP = info[0]
		if currentSNP%5000 == 0:
			print('SNP {0} of {1} SNPS.  Working...'.format(currentSNP, lastSNP))
		infoList = [] #col0=SNPnum, col1=SNPallele, col2 = SNPpos, col3 = GenomeID, col4 = AccNum, col5=geneType, col6 = codon, col7=aa, col8=peptide, col9 = product
		InAnnotatedRegion = 'F'
		InAnnotatedGenome = 'F'
		
			
		#********* end block 1 ********
		#********* start block 2 ******
	SNPnum = info[0]
	SNPallele = info[2]
	GenomeID = info[3]
	SNPpos = info[4]	
	codonInfo = {} #key is SNPallele value is annonymous list [codon, aa]
	if SNPpos != 'x':		
		locusSeq = info[1]
		
				
		strand = info[5]
		AccNum = info[6]
		if AccNum != 'NoAccNum':
			InAnnotatedGenome = 'T'
			InAnnotatedRegion = 'F'
			for i in range(len(RefGenomesDict[AccNum])):
				if SNPpos >= RefGenomesDict[AccNum][i][2] and SNPpos <=RefGenomesDict[AccNum][i][3]:
					geneType = RefGenomesDict[AccNum][i][0]
					complement = RefGenomesDict[AccNum][i][1]
					start = RefGenomesDict[AccNum][i][2]
					end = RefGenomesDict[AccNum][i][3]
					product = RefGenomesDict[AccNum][i][4]
					if geneType == 'CDS':
						InAnnotatedRegion = 'T'
						(peptide, aa, codon, readingFrame,RC) = findAA(locusSeq, SNPallele, SNPpos, start, end, strand, complement)
						infoList.append([SNPnum, SNPallele, SNPpos, GenomeID, AccNum, geneType, codon, aa, peptide, product])
						if SNPallele not in codonInfo.keys():
							codonInfo[SNPallele] = [codon, aa, geneType, peptide, product]						
						#print("SNP {0} is bp {7} in a CDS that is bp {1} through {2} in {3}\n It encodes {4} in the peptide {5} in {6}. ".format(SNPnum, start, end, GenomeID, aa, peptide, product, SNPpos))
						break
					elif geneType == 'RNA'or geneType == 'tRNA':
						InAnnotatedRegion = 'T'
						infoList.append([SNPnum, SNPallele, SNPpos, GenomeID, AccNum, geneType, '?', '?', '?', product,RC])
						if SNPallele not in codonInfo.keys():
							codonInfo[SNPallele] = ['?', '?', geneType, '?', product,RC]
				else:
					if 'InAnnotatedRegion' in locals() and InAnnotatedRegion == 'T':
						pass
					else:
						InAnnotatedRegion = 'F'
			if InAnnotatedRegion == 'F':
				infoList.append([SNPnum, SNPallele, SNPpos, GenomeID, '?', '?', '?', '?', '?', '?', '?'])
	else:  #if in an unannotated genome, i.e. SNPpos is 'x'
		if 'InAnnotatedGenome' in locals() and InAnnotatedGenome == 'T':
			pass
		else:
			InAnnotatedGenome = 'F'
			infoList.append([SNPnum, SNPallele, SNPpos, GenomeID, '?', '?', '?', '?', '?', '?'])
	#********* end block 2 ******
OUTFILE.close()
SUMMARY.close()
endTime = time.time()
usedTime = endTime - startTime
print("\nUsed {0:.1f} seconds to parse annotations.".format(usedTime))
maxUsage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print("Max memory used to parse annotations was {0} kilobytes.".format(maxUsage/1024))

