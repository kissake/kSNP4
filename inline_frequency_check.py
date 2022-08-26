#!/usr/bin/python3

# This script is intended to be used in a pipeline to calculate statistics on a file as
# it is being written to disk.  You might use it to process the significant output of
# a 'jellyfish dump' run without having to re-read the data from disk.  Example:
#
# jellyfish dump -c Jelly.fsplit0 | inline_frequency_check freq.fsplit0 >> unsortedkmers.fsplit0
#
# The statistics output will be in the file named: freq.fsplit0.
#
#
# Note that this may also be imported as a library, and be similarly used with a pipe.


import statistics
import sys


def getMinKmerCoverage(inFile, outFile):

    '''
    The common case for this function is a lot of kmers with frequency 1
    and a small number (under 1%) with a larger frequency, typically under
    100 even in those cases.  Stats functions on lists are _slow_ for large
    lists, so there are some sortcuts we take math-wise to get a faster
    answer without compromising accuracy.

    Without this optimization, calculating the statistics takes about as
    long as passing the file.  With the optimization, the stats do not
    take an easily measureable amount of time.
    '''
    
    countsList = []
    space = " "
    # Some simple math to speed mean and median calculation in typical cases.
    # If you subtract one from each element in a series, the average of that
    # series will be one less than the average of the original series.
    runningSumMinusOne = 0
    # This is a count of the frequency of occurrences of the number 1 in
    # the data.
    onesCount = 0


    for line in inFile:
        
        # Later, process lines in batches of N where N is big, to see if
        # we can get a performance boost.  readlines() takes a "hint" argument
        # that lets us specify a maximum size of batches of lines to read....
        
        outFile.write(line)
        merCount = int(line.split(space)[1])
        # No sense appending a '1' to this list... TODO FIXME XXX - JN
        countsList.append(merCount) #Split line on spaces, get second element.
        if merCount == 1:
            onesCount = onesCount+1
            # no need to add anything to runningSumMinusOne; 1-1 = 0
        else:
            runningSumMinusOne = runningSumMinusOne + merCount - 1


    # Fastpath for the common case: Typically, far more than 50% of the values are 1,
    # so no need to parse the entire array to come up with this result.
    if onesCount > (len(countsList) / 2):
        # If more than half of the values are a particular value, the median is that
        # value by definition.
        median = 1
        # Definition of the mean: The sum divided by the count.
        # Note that python3 integers _cannot overflow_ (you will get OOM instead)
        mean = ( runningSumMinusOne / len(countsList) ) + 1
    else:
        # THIS PATH IS SLOW.  At least two O(N) operations, possibly O(NlogN)
        # Math on countslist
        mean = statistics.mean(countsList)

        # NOTE: This operation is _SLOW_ (at least O(N).  Would be better to find a way to compute
        # the running median?  Keep a sorted array?
        median = statistics.median(countsList)


    # Halfway between median and mean)
    return int(statistics.mean([mean, median]))

            
if __name__ == "__main__":

    # Variables and constants:
    newline = "\n"

    outputFilename = sys.argv[1]
    outputFile = open(outputFilename, 'w')
    outputStat = getMinKmerCoverage(sys.stdin, outputFile)
    outputFile.close()

    sys.stdout.write(str(outputStat) + newline)
