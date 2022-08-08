#!/usr/bin/python3

import sys

# Previous versions of this program accepted a value for k on the commandline and
# worked out an appropriate value for prefixLen.


# prefixLen of 4 means 256 different files.  Some OSes have a limit on number of open files, e.g. 1k
# If this value needs to be changed, consider testing for # of entries in buckets (len(buckets)), and
# when it gets too big, add logic for closing the file for the first bucket (buckets[buckets.keys()[1]].close()) and then deleting it
# del buckets[buckets.keys()[1]]
# At this point, benefit may be gained from a sorted input (but not necessary; just means closing and
# opening files more than otherwise, so a performance hit)

prefixLen = 4  
outputBuffer = 100000


def subset4mers(inputfile, prefixlen):
    # Holds the file handle for each file we have opened, one per bucket.
    buckets = {}

    in = open(inputfile, "r")

    line = in.readline()
    while line:
        bucketName = line[:prefixlen]

        # Open a file if we don't have one open for this prefix.  Otherwise
        # return the file we already opened for this prefix.
        bucketFile = buckets.setdefault(bucketName, open(inputFile + "." + bucketName, "a"))
        bucketFile.write(line)
  
        line = in.readline()
    return buckets
        
if __name__ == "__main__":

    inputFile = sys.argv[0]
    if prefixLen <= 4:
        subset4mers(inputFile, prefixLen)

    else:
        print("Prefix lengths greater than 4 not supported (yet)")
        # suggestion for implementation is to modify subset4mers to take an additional
        # argument indicating how many chars to skip to get to the chars to bucket, then
        # run iteratively until you get to the desired length.  Requires more storage and
        # lifting, may not be worth avoiding sorted input at that point?
        sys.exit(1) 
    
