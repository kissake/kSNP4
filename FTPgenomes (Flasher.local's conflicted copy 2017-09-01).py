#!/usr/bin/env python

"""Description: 
infile is list of genomes, one genome per line in format ID type path
where ID is genome name that may contain spaces, type is complete, scaffold or contig,
and path is path to directory containing the .fna file
Usage: FTPgenomes infileName

"""
#********************* MODULES**********************
from __future__ import division
from __future__ import print_function
import sys
import os
import time
import ftputil
import gzip #permits extracting the .gz files from draft genomes
#******************* FUNCTION DEFINITIONS *******************
def getFile(ID, path):
	#get the name of the gz file to be retreived
	temp = path.split('/')
	gzFile = temp[-1]+'_genomic.fna.gz'  #gzFile is the gzipped file to be downloaded
	
	#get the directory from which to retrieve the file
	temp = path.split('/')
	del temp[0:3]
	dir = '/'.join(temp)
	print(dir)

	#connect to the source database
	host = ftputil.FTPHost('ftp.ncbi.nih.gov', 'anonymous', 'password')

	#open the directory given in path
	host.chdir(dir)

	#download the file
	print("Downloading " + gzFile + ' ...')
	host.download(gzFile,gzFile) 
	
	#extract the .gz archive and save it as ID.fna
	archive = gzip.open(gzFile)
	outfileName = ID+'.fna'
	OUTFILE = open(outfileName, 'wb')
	OUTFILE.write(archive.read())
	archive.close()
	os.remove(gzFile)

	host.close()
#********************** Main *************************
startTime = time.time()

infile = sys.argv[1]
INFILE = open(infile, "rU")
for line in INFILE:
	line = line.rstrip()
	temp = line.split('\t')
	ID = temp[0]
	path =temp[2]
	getFile(ID, path)


endTime = time.time()
elapsedTime = endTime-startTime

print("Used {0} seconds to download these files.\n".format(elapsedTime))
