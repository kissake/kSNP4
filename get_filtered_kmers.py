#!/usr/bin/python3

'''
Expected arguments:
 - k (kmer length)
 - fastalist (input file with 
'''

import argparse
import re
import subprocess
import os
import sys # Replace with logging!
import logging as logging
import hashlib
import ksnpCache as kcache
import shutil



tab = '\t'
newline = '\n'

##################################
##################################
#
# FUNCTIONS
#
##################################
##################################


def parseGetFilteredKmerArgs( override = None ):
    parser = argparse.ArgumentParser(description='output a fileName2genomeName file and corresponding data files')
    parser.add_argument('fastalist', help='The list of fasta files used as input to kSNP')
    parser.add_argument('mergeFastaReads', help='Command to merge contigs into a single fasta file')
    parser.add_argument('k', help='The kmer length for this run')
    parser.add_argument('jellyfish', help='The path to the jellyfish program to use for this run')
    parser.add_argument('hashSize', help='The hash size to pass to the jellyfish program')
    parser.add_argument('numCpus', help='The number of CPUs to use for the jellyfish run.')
    parser.add_argument('fileToGenome', help='The filename of the fileName2genomeName output')
    parser.add_argument('freqCheck', help='The filename of the program to calculate count frequencies inline')
    parser.add_argument('outputPrefix', help='The prefix to use for output files, if any', default = 'fsplit')
    parser.add_argument('kmerSuffix', help='The suffix to use for files cotaining kmer frequency counts', default = 'kmer')
    parser.add_argument('jellyfishSuffix', help='The suffix to use for files containing jellyfish database data', default= 'Jelly')
    parser.add_argument('filterSuffix', help='The suffix to use for filtered kmer frequency counts', default= 'filtered')
    kcache.addCacheOptions(parser)
    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')
    parser.add_argument('--info', action='store_true', help='Output progress updates to the screen using STDERR')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)


def mergeFastaReads(inputFile, inputFileName, outputFile):
    '''
    BROKEN
    '''
    
    filler = 'N'
    headerMatch = re.compile(r'^>.*$') # To be removed
    firstHeaderRE = r'^>.*' + newline
    subsequentHeadersRE = newline+'>.*'+newline
    headerMatch = re.compile( firstHeaderRE + r'|' + subsequentHeadersRE)
    # Match two or more Ns (if there is only one N, no substitution to do.
    multipleNsMatch = re.compile(r'NN+')
    # A str.translate table to remove newlines from a string. (faster than re)
    removeNewlines = str.maketrans(newline,filler)

    

    #Output
    outputFile.write(">" + inputFileName + " merged" + newline)

    # No idea how efficient this is.  It could blow up memory?
    outputFile.write(
        multipleNsMatch.sub(
            filler, headerMatch.sub( # Remove multiple Ns.
                '',inputFile.read() # Remove header lines
                ).translate(removeNewlines)
            )
        )
    outputFile.write(newline)
    


def parseFastaInputFiles(fileToGenomeFilename,fastalistFilename, inputFilenamePrefix, jellyfishFilenameSuffix, kmerFilenameSuffix, k):

    '''
    Given input filename, and information about the expected filenames, create the
    fileName2genomeName file and set the names for the temporary (in-processing)
    per-gene files.  Return a list of input files and important details about them.
    '''
    
    inputLineCounter = 0
    # List of lists, each element is the number, genome, and filename for an input file)
    inputFiles = []

    logging.debug("Parsing input file %s to generate %s and calculate filenames", fastalistFilename, fileToGenomeFilename)
    fileToGenomeFile = open(fileToGenomeFilename, 'w')

    # Assume this is manageably small.  Typically 100 entries would be considered huge.
    inputLines = open(fastalistFilename, 'r').readlines()
    for line in inputLines:
        # For each line in the input list of fasta files...
        
        (filename, genome) = line.split(tab)
        countStr = str(inputLineCounter)

        logging.debug("Calculating checksum for %s",filename)
        # Calculate the checksum to determine the filename for this data in the cache.
        thisFile = open(filename,'rb')
        fileHash = hashlib.sha256()
        fileHash.update(thisFile.read())
        checksum = str(fileHash.hexdigest())
        thisFile.close()
        logging.debug('Determined checksum for %s: %s', filename, checksum)

        thisRecord = {
            'count': countStr,
            'genome': genome.strip(),
            'inputFile': filename,
            'tempFile': inputFilenamePrefix + countStr,
            'jellyfishDB': inputFilenamePrefix + countStr + '.' + options.jellyfishSuffix,
            'kmers': inputFilenamePrefix + countStr + '.' + options.kmerSuffix,
            'kmersFiltered': options.kmerSuffix + '.' + inputFilenamePrefix + countStr, # Set by expectations of the rest of the parent script
            'checksum': checksum,
            'cacheId': cacheFileName(checksum, k),
        }

        # I would prefer the kmersFiltered name be as defined below to more accurately represent origin / processing.
        # 'kmersFiltered': inputFilenamePrefix + countStr + '.' + options.kmerSuffix + '.' + options.filterSuffix,
        # Even better would be:
        # 'kmersFiltered': '.'.join([genome.strip(), checksum[-10:], options.kmerSuffix, options.filterSuffix]) 
        
        fileToGenomeFile.write(thisRecord['tempFile'] + tab + thisRecord['genome'] + newline)


        # add this entry to our list of files to process
        inputFiles.append(thisRecord)
        inputLineCounter = inputLineCounter + 1
        
    fileToGenomeFile.close()

    logging.info("Number of input sequences: %s", len(inputFiles))

    return inputFiles

def getCachedFilteredKmers(fileData, ):
    # If there is cached data for these files, copy the cached data into the final location
    # and return a list of the files for which we did NOT find cached data (or an empty list
    # if we got it all)
    return fileData


def convertInputFiles(inputFileList, mergeFastaReadsFilename):
    for record in inputFileList:
        print(tab.join([record['count'], record['genome'], record['inputFile']]))
        outputFileHandle = open(record['tempFile'], 'wb')
        
        subprocess.run([mergeFastaReadsFilename, record['inputFile']], stdout=outputFileHandle)
        outputFileHandle.close()

def processJellyfishDumps(inputFiles, jellyfishFilename, freqcheckFilename, ):
    # Now we grab frequency data from these databases.
    # Note that we run 'dump' in the background, so if there are a v. large number of input files,
    # this could get exciting. (but probably not; we'll just go a little slower due to resource
    # contention, but still be much faster because we're letting the OS manage it, rather than one-
    # by-one)
    
    allPendingProcesses = []
    for file in inputFiles:
        
        nextFile = file['jellyfishDB']
        dbFileCount = 0
        processes = []

        logging.debug("Processing (dump) %s", file['jellyfishDB'])
        
        while os.access(nextFile, mode=os.F_OK):

            if dbFileCount > 0:
                # This will fail to detect the case where there is fsplitN.Jelly and fsplitN.Jelly_0.
                logging.warning("WARNING: The case of multiple jellyfish dump outputs is not handled at this time.")

            
            jellyfishDump = subprocess.Popen([
                jellyfishFilename, 'dump', '-c', nextFile], stdout = subprocess.PIPE)
            
            freqCheck = subprocess.Popen([
                freqcheckFilename, file['kmers'], ], stdin=jellyfishDump.stdout, stdout=subprocess.PIPE)

            # Permit SIGPIPE to reach jellyfish dump process by ensuring that we don't hold an open
            # file handle to the pipe in question.
            jellyfishDump.stdout.close() 

            processes.append(freqCheck)
            allPendingProcesses.extend([jellyfishDump, freqCheck])

            # Join STDOUT of jellyfishDump to STDIN of freqCheck, and capture STDOUT of freqCheck.

            # The intent is that nextFile always be 1 behind dbFileCount so we can hit the un-suffixed file
            # first and catch up with _0.
            nextFile = file['jellyfishDB']+ '_' + str(dbFileCount)
            dbFileCount = dbFileCount+1

        file['freqProcesses'] = processes

    # Wait for all of our dumpers and frequency calculators to finish
    logging.debug("Waiting for dumpers and frequency checkers to finish")

    # Busy wait for processes starting at 0.
    while len(allPendingProcesses) > 0:
        # This polling is ONLY okay because we know the output from inline_frequency_check is
        # v. tiny (integer representation smaller than size of the source file), and won't fill
        # the pipeline.  OTHERWISE, we would need to find a way to grab the data so the process
        # doesn't wait to output while we wait for it to quit (deadlock).
        pollResult = allPendingProcesses[0].poll()
        if pollResult is not None:
            allPendingProcesses.pop(0)
            logging.debug('a process exited with exit code %s, %s processes remaining to check', pollResult, len(allPendingProcesses))
            if pollResult != 0:
                raise RuntimeError ("Subprocess exited with nonzero exit code ERROR!")
    logging.debug("Dumpers and statistics done")

    # Grab frequency data from per-process output.
    for file in inputFiles:

        logging.debug('Processing frequency data for %s', file['tempFile'])
        # Grab the frequency cutoff from the first subprocess, ignoring stderr.
        processes = file['freqProcesses']
        file['minFreq'] = int(processes.pop().communicate()[0])

        # Then clean up any remaining processes for this file (shouldn't be any)
        for process in processes:
            process.communicate()
            
        print(file['tempFile'] + tab + str(file['minFreq']))

    # The only data relevant in this data structure is 'minFreq', but this is easy.
    return inputFiles
        

def generateJellyfishDumps(inputFileList, jellyfishFilename, k, hashSize, numCPUs):
    '''
    Process all of the input files to get lists of kmers of all frequencies.
    '''
    
    # Jellyfish is pretty resource intensive; run these processes one-at-a-time.
    # Further investigation suggests that for moderate sized inputs (~200M), peak resource
    # usage is only for the first 10 seconds or so... So can we run multiple at the same time?
    
    for file in inputFileList:

        logging.debug("Processing (jellyfish) %s", file['jellyfishDB'])

        # Note that this output is e.g. file0.Jelly, or file0.Jelly_<n>
        subprocess.run([
            jellyfishFilename,
            'count', '-C', '-o', file['jellyfishDB'], '-m', k, '-s', hashSize, '-t', numCPUs, file['tempFile'] ])

    # Determine if we need to do any merging by looking for the clearest sign of a second dump file: <output>_1.
    logging.debug("Checking for multiple jellyfish output files; this is new; previously we just processed each individual output file.")
    for file in inputFileList:
        if os.access(file['jellyfishDB'] + '_1', mode=os.F_OK):
            # At least two dump files exist, they need to be merged
            logging.warning("WARNING: If this file has a min_kmer_frequency greater than 1, there is the chance to miss some SNPs.")
            logging.warning("WARNING: It will probably work okay if you comment out this following 'raise', but see comments in source.")

            # NOTE NOTE NOTE:  Ignoring this is technically WRONG in the case where there is more than one
            # output file.  Failure mode is:
            #  - raw reads (so our min_kmer_freq is higher than 1)
            #  - kmer spread between multiple files and not notably well represented.
            # Correct behavior is to jellyfish merge.  FIXME TODO XXX - JN
            # I think the right thing to do for now is to output a warning to the user.
            
            raise ValueError ("Jellyfish output overflowed.")

def cacheFileName(checksum, k):
    return (checksum + "." + str(k))
        
def cacheFilteredKmers(inputFiles, k, cacheDir):
    if cacheDir == '':
        # Caching disabled.
        logging.info('Caching disabled due to --cachedir not set / set to the empty string')
    else:
        for file in inputFiles:
            id = file['cacheId']
            logging.debug('Cached filtered kmers from %s to %s.', file['kmersFiltered'] , kcache.cacheFileName(id, cacheDir, kcache.filteredKmersData))

            shutil.copy(file['kmersFiltered'], kcache.cacheFileName(id, cacheDir, kcache.filteredKmersData))
        # TODO FIXME DEBUG XXX We should shrink the cache to size here. - JN
        
        

##################################
##################################
#
# MAIN
#
##################################
##################################

    

if __name__ == "__main__":
    uncachedFiles = []
    
    options = parseGetFilteredKmerArgs()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    elif options.info:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARNING)

    logging.debug('Commandline options, as parsed: %s', str(options))

    logging.info('Parsing input fasta file')
    inputFiles = parseFastaInputFiles(options.fileToGenome,options.fastalist, options.outputPrefix, options.jellyfishSuffix, options.kmerSuffix, options.k)

    if options.cachedir == '':
        # Caching disabled
        logging.info('Caching disabled, processing all %s input files.', len(inputFiles))
        uncachedFiles = inputFiles
    else:
        # Cache specified.
        logging.info('Checking cache for %s input files', len(inputFiles))
        for file in inputFiles:
            cacheFilename = kcache.cacheFileName(file['cacheId'], options.cachedir, kcache.filteredKmersData)
            if os.access(cacheFilename,mode=os.R_OK):
                # if the data is cached, copy it into place (why do math; it has already been done)
                logging.info('%s found in cache, restoring cached data and skipping processing.', file['inputFile'])
                shutil.copy(cacheFilename, file['kmersFiltered'])
            else:
                # We didn't find this data in the cache.
                uncachedFiles.append(file)

    if len(uncachedFiles) > 0:
        # Convert input files into format appropriate for kSNP processing.
        logging.info('Converting %s input files without cached data into canonical kSNP input', len(uncachedFiles))
        convertInputFiles(uncachedFiles, options.mergeFastaReads)
        
        # Run Jellyfish on the input files to produce kmers list files
        logging.info("Running jellyfish to find kmers for %s input files", len(uncachedFiles))
        
        generateJellyfishDumps(uncachedFiles, options.jellyfish, options.k, options.hashSize, options.numCpus)
        dumpFiles = processJellyfishDumps(uncachedFiles, options.jellyfish, options.freqCheck)

        logging.debug("Filtering out minimum frequency kmers from %s files", len(dumpFiles))

        for file in dumpFiles:
            logging.debug("File: %s", file)
            
            if file['minFreq'] != 1:  #Only if there is any filtering to do:
                logging.debug("Running awk to filter minimum frequency kmers (awk-ward) for %s", file['tempFile'])
                kmersFiltered = open(file['kmersFiltered'], 'wb')
                subprocess.run(['awk', '-v', 'm='+str(file['minFreq']), '$2>=m {print}', file['kmers']], stdout=kmersFiltered)
                kmersFiltered.close()
            else:
                logging.debug("No filtering to do, renaming kmers file to filtered kmers filename for %s", file['tempFile'])
                os.rename(file['kmers'],file['kmersFiltered'])

        # Cache results (if possible; if cachedir set to empty string, this is a noop.
        logging.debug("Caching results.  Files to cache: %s", dumpFiles)
        cacheFilteredKmers(dumpFiles, options.k, options.cachedir)

    else: # We didn't find any uncached files
        logging.info('Cached data found for all %s input files, using cached data only.', len(inputFiles))
        
    exit(0)