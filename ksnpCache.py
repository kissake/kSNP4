#!/usr/bin/python3

import logging as logging

# Import standard Python libraries for argument parsing, interacting with the
# filesystem, and environment and time / date processing.
import argparse as argparse
import os as os
import datetime as datetime


# Static

# Default values absent other sources of defaults (i.e. initial value for placing into
# a config file, or value to use if caching is enabled, but these values not specified.
defaultMaxCacheBytes=10000000000
defaultMaxCacheDays=90
defaultCacheDir=os.path.join(os.environ.get('HOME',''),'kSNP')

annotationData = 'gb'
genomeData = 'fasta'
filteredKmersData = 'kmers.filtered'
cacheLocation = {
    annotationData: 'GbkCache',
    genomeData: 'GenomeCache',
    filteredKmersData: 'FilteredKmersCache',
}

def parseCommandline(override=None):

    example='''To get the current size of the cache:

%(prog)s --stats'''
    
    parser = argparse.ArgumentParser(description='This is a module for getting information about and manipulating the kSNP cache.',
                                     epilog=example, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--stats', action='store_true', help="Show the current cache statistics, including size.")
    parser.add_argument('--shrink', action='store_true', help="Shrink the cache to the current maximum size.")
    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')

    # Add standard options for cache-enabled programs
    addCacheOptions(parser)
    
    parser.add_argument('--annotation', dest='responses', action='append_const', const=annotationData,
                        help='Output information about annotation data.')
    parser.add_argument('--genome', dest='responses', action='append_const', const=genomeData,
                        help='Output information about genome data.')
    parser.add_argument('--filteredKmers', dest='responses', action='append_const', const=filteredKmersData,
                        help='Output information about cached filtered kmers.')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)


def addCacheOptions(parser):
    '''
    A utility for programs using this library to add options to their argument parser that will permit the user
    to set values for and take advantage of defaults for various caching related settings.
    '''
    # Longer term, this will help us set these using environment variables or config files without having to modify dependent programs.
    
    # The default will eventually be configurable with a config file.
    parser.add_argument('--cachedir', default=defaultCacheDir,
                        help='Override the directory to use for caching data')
    # The default will eventually be configurable with a config file.
    parser.add_argument('--maxcachedays', type=int, default=defaultMaxCacheDays, help='Override the maximum age for cached data in days')
    parser.add_argument('--maxcachebytes', type=int, default=defaultMaxCacheBytes, help='Override the maximum size for the cache in bytes')


def cacheFileName(accession, cacheDir, dataType=annotationData):
    # Combine the base cache area plus the specific data type being cached with the accession id.

    if cacheDir == '':
        # Caching disabled.  cacheFileName isn't meaningful.
        raise ValueError("Caching disabled, cannot generate filename.")
        return ''
        
    else: # Caching enabled, return answer.
        cacheFullPath = os.path.join(cacheDir, cacheLocation[dataType])
        
        # Test to see if the cache directory exists.
        if not os.access(cacheFullPath, mode=os.F_OK):
            # If the cache directory doesn't exist, create it.
            os.makedirs(cacheFullPath)
            
        return os.path.join(cacheFullPath, str(accession))


def cacheSize(cacheDir, dataType=None):
    # Default is to get the total size of the cache, in bytes.
    size = 0
    sizeDir = cacheDir

    if cacheDir == '':
        # Caching disabled.  Cache size is meaningless.
        raise ValueError("Caching disabled, cannot generate filename.")
        return ''
        
    else: # Caching enabled, return answer.
        if dataType is not None:
            sizeDir = os.path.join(cacheDir, cacheLocation[dataType])
            logging.debug("Showing size data for %s in %s", dataType, sizeDir)
            
            
        for root, dirs, files in os.walk(sizeDir):
            fileSizes = [os.path.getsize(os.path.join(root, name)) for name in files]
            size = size + sum(fileSizes)
                
        return size

def shrinkCache(cacheDir, newSize):
    # Shrink the size of the cache to a size less than the new size by deleting
    # files that are the most stale.

    startingSize = 0
    allFiles = []

    if cacheDir == '':
        # Caching disabled.  cacheFileName should not be called.
        raise ValueError("Caching disabled, cannot generate filename.")
        return ''
        
    else: # Caching enabled, return answer.

        logging.debug("Shrinking cache to %s bytes", newSize)
        
        for root, dirs, files in os.walk(cacheDir):
            for name in files:
                fullName = os.path.join(root, name)
                data = os.stat(fullName)
                fileDbEntry = (fullName, data.st_atime, data.st_mtime, data.st_ctime, data.st_size)
                startingSize = startingSize + fileDbEntry[4]
                allFiles.append(fileDbEntry)
                
        if startingSize < newSize:
            # We already achieved the goal, no need to delete files.
            logging.debug("Target size (%s) larger than current cache size(%s).", newSize, startingSize)
            return ( 0, startingSize )
        else:
            # Order the files in the cache by oldest (least recently used) to newest, then
            # progress through the list until we have achieved the target size.
            sizeChangeRequired = startingSize - newSize
            cumulativeSize = 0
            filesToDelete = []
            allFiles.sort(key=lambda entry: entry[1]) # Sort on atime first.
            raise NotImplementedError("This path not tested.  May delete newest cached files.")
            logging.debug("Oldest 10 files: %s", allFiles[:10])
            logging.debug("Newest 10 files: %s", allFiles[-10:])
            while cumulativeSize < sizeChangeRequired:
                filesToDelete.append(allFiles.pop())
                cumulativeSize = cumulativeSize + filesToDelete[-1][4]

            # Delete the files required to reach the required size.
            logging.debug("Deleting %s files", len(filesToDelete))
            for file in filesToDelete:
                os.remove(file[0])

            return (len(filesToDelete), startingSize - cumulativeSize)
        
            


##############################################################################
##############################################################################
#
# Main / cache management functions
#
##############################################################################
##############################################################################

if __name__ == "__main__":

    options = parseCommandline()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s')

    logging.debug('Commandline options, as parsed: %s', str(options))

    if options.cachedir == '':
        print("Caching disabled (cachedir set to empty string).  Nothing to report.")
        exit(0)
        
    if options.stats:
        if options.responses is not None:
            for response in options.responses:
                print("Size of data stored in %s portion of cache: %s" % (cacheLocation[response], cacheSize(options.cachedir, dataType=response)))
               
        print("Total cache size: %s" % ( cacheSize(options.cachedir) ) )
        
    if options.shrink:
        (numDeleted, newSize) = shrinkCache(options.cachedir, options.maxcachebytes)
        print("Deleted %s files, shrinking cache to %s (less than or equal to than maximum: %s)" % (numDeleted, newSize, options.maxcachebytes))

    if not ( options.shrink or options.stats ):
        # One of these is required
        parseCommandline(['--help'])

