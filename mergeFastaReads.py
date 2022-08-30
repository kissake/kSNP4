#!/usr/bin/python3

'''
Expected arguments:
 - k (kmer length)
 - fastalist (input file with 
'''

import re
import sys # Replace with logging!
import logging as logging
import hashlib as hashlib



def mergeFastaReadsWithChecksum(inputFile, inputName, outputFile):
    # This takes more than 1/3 longer than without checksum (doesn't take
    # advantage of multi-core), but it is still faster than the perl, and
    # less than 50% worse than the checksum alone.
    # 
    # Checksum alone, commandline tool: real	0m2.963s
    # mergeFastaReads alone:            real	0m3.235s
    # This function (both):             real	0m4.425s
    # Legacy mergeFastaReads:           real	0m7.754s
    

    # This is a moderately productive multiple of 4096, which is a typical
    # block size.
    buffersize = 32768
    firstHeader = True
    checksum = hashlib.sha256()

    # Write the fasta header for the output file.
    outputFile.write(">" + inputName + " merged\n")

    # Fill the buffer
    data = inputFile.read(buffersize)
    checksum.update(data.encode('utf-8'))
    bytesRead = len(data)

    # Test for end-of-file (no data read / remaining)
    while bytesRead > 0:
        (remainingData, output) = processBuffer(data)
        if firstHeader and output != '' and output[0] == 'N':
            firstHeader = False
            # Omit the initial N
            outputFile.write(output[1:])
        else:
            outputFile.write(output)
        newData = inputFile.read(buffersize)
        checksum.update(newData.encode('utf-8'))
        bytesRead = len(newData)
        data = remainingData + newData


    # Process the remaining data, making sure to drop any lingering header data
    (remainingData, output) = processBuffer(data, last=True)
    # Write _all_ of the output data including any trailing 'N's, and a terminating
    # newline (we got rid of almost all of the others)
    outputFile.write(output + remainingData + '\n')

    return checksum.hexdigest()


def mergeFastaReads(inputFile, inputName, outputFile):

    # This is a moderately productive multiple of 4096, which is a typical
    # block size.
    buffersize = 32768
    firstHeader = True

    # Write the fasta header for the output file.
    outputFile.write(">" + inputName + " merged\n")

    # Fill the buffer
    data = inputFile.read(buffersize)
    bytesRead = len(data)

    # Test for end-of-file (no data read / remaining)
    while bytesRead > 0:
        (remainingData, output) = processBuffer(data)
        if firstHeader and output != '' and output[0] == 'N':
            firstHeader = False
            # Omit the initial N
            outputFile.write(output[1:])
        else:
            outputFile.write(output)
        newData = inputFile.read(buffersize)
        bytesRead = len(newData)
        data = remainingData + newData


    # Process the remaining data, making sure to drop any lingering header data
    (remainingData, output) = processBuffer(data, last=True)
    # Write _all_ of the output data including any trailing 'N's, and a terminating
    # newline (we got rid of almost all of the others)
    outputFile.write(output + remainingData + '\n')



def processBuffer(data, last=False):
    '''
    Take data, and return it in two sections:
     - Output: Any output that could have been derived from the data, with the corresponding
         data consumed (not returned in the queue)
     - Queue: Any unconsumed data.

    This function is responsible for turning header lines into Ns.
    '''
    
    headerChar = '>'
    newline = '\n'
    N = 'N'
    queue = ''
    output = ''

    dataLines = data.split(newline)
    if last:
        # Accommodates the case where the last line is a header line and does not
        # have a terminating newline.  This is bad input in two different ways.
        dataLines.append('')

    for nextLine in dataLines[:-1]:
        # If we make it into the loop, the current line has a corresponding newline
        # and can be processed fully.
        if nextLine == '' or nextLine[0] != headerChar:
            queue = queue + nextLine
        else:
            queue = queue + N

    nextLine = dataLines[-1]
            
    if nextLine == '' or nextLine[0] != headerChar:
        # This isn't a header line, we can just process this data; we will catch up
        # still in the middle of this line.
        
        return processQueue(queue+nextLine)
    else:
        # This last line is a header line, but we don't have a terminating newline.
        # We may have a bunch of data to process, and we'll have to prepend any
        # pending data in front of nextline, so to represent reality accurately,
        # we will put a newline between pending data and this start-of-header-
        # line.

        # Just process the data that is pending, not including this line.
        (queue, output) = processQueue(queue)

        # any remaining queued data is prepended to this line, with a newline
        # between:
        return (queue + newline + nextLine, output)





def processQueue(queueData):
    '''
    Take data without header information and squash repeating N's to prepare for
    output.
    '''

    newline = '\n'
    repeatNs = r'NNN*'
    
    # We know there are no headers in this business, just filter
    # Remove newlines
    cleaned = queueData.translate(str.maketrans('','',newline))
    squashed = re.sub(repeatNs, 'N', cleaned)
    if squashed != '' and squashed[-1] == 'N':
        return ('N', squashed[:-1])
    else:
        return ('', squashed)

    
        
if __name__ == "__main__":
    inputFile = open(sys.argv[1], 'r')

    # mergeFastaReads(inputFile, sys.argv[1], sys.stdout)
    sys.stderr.write(mergeFastaReadsWithChecksum(inputFile, sys.argv[1], sys.stdout))

    
