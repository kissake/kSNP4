#!/usr/bin/env python

"""
Description: Writes a kSNP4 infile for the files in the target directory. Must be run from 
within the directory ENCLOSING the target directory

Usage: MakeKSNP4infile -indir Target directory -outfile name of outfile [S] where
S is an optional argument for the semi-automatic mode.  The mode is A for automatic by default
"""
#********************* MODULES**********************
import sys
import os
#******************* FUNCTION DEFINITIONS *******************



#********************** Main *************************
#variables
mode = 'A' #A is the automtic mode, the alternative is S
indir = ''
outfile = ''
home = os.getcwd()
thisPath = ''
genomeID = ''
Files = [] #list of files in target directory
fileList = [] #col0 is path col1 is file name

#get command line arguments
for i in range(1, len(sys.argv)):
	if sys.argv[i] == '-indir':
		indir = sys.argv[i+1]
	if sys.argv[i] == '-outfile':
		outfile = sys.argv[i+1]
		if sys.argv[i] == 'S':
			mode = 'S'

os.chdir(indir) #make the target directory the cwd
thePath = os.getcwd()
#print(thePath)
#print('After changing to the input directory the path is now {0}\n'.format(thePath))



#get the list of files
Files = os.listdir('.')
for file in Files:
	if not file.startswith('.'):
		#sys.stdout.write("{0}\n".format(file))
		thePath = os.path.abspath(file)
		#OUTFILE.write('{0}\t{1}\n'.format(thePath, file))
		fileList.append([thePath, file])
#print(Files[0])

os.chdir(home)
OUTFILE = open(outfile, 'w')
for i in range(len(fileList)):
    OUTFILE.write('{0}\t{1}\n'.format(fileList[i][0], fileList[i][1]))


OUTFILE.close()
