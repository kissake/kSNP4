#!/usr/bin/env python3

'''
usage: getGenomeDataByAccession.py [-h] [--debug] [--cachedir CACHEDIR] [--maxcachedays MAXCACHEDAYS] [--annotation] [--genome] accessionID [accessionID ...]

Return genome data for given accession IDs

positional arguments:
  accessionID           Accession IDs to return.

optional arguments:
  -h, --help            show this help message and exit
  --debug               Output diagnostic data to the screen using STDERR
  --cachedir CACHEDIR   Override the directory to use for caching data
  --maxcachedays MAXCACHEDAYS
                        Override the maximum age for cached data
  --annotation          Output annotation data for the accession ID if available.
  --genome              Output genome data for the accession ID if available.

To get the genome (FASTA format) data for accession ID NZ_CP075697.1 to your screen:

getGenomeDataByAccession.py --genome NZ_CP075697.1

'''


# Import the Entrez database API from BioPython
import Bio.Entrez as Entrez

import ksnpCache as kcache

# Import standard Python libraries for argument parsing, interacting with the
# filesystem, and environment and time / date processing.
import argparse as argparse
import os as os
import datetime as datetime

# Permit use of stdout
import sys

# Help with debug data:
import logging as logging

# Import tools for parsing kSNP configuration data.
# TODO - JN
# import ksnpConfig


####################################
####################################
###
###
### Functions
###
###
####################################
####################################


# The override option to this function can be used to cause the parsing of a provided list
# instead of the programs passed options, for testing.
def parseCommandline(override=None):

    example='''To get the genome (FASTA format) data for accession ID NZ_CP075697.1 to your screen:

%(prog)s --genome NZ_CP075697.1
'''
    
    parser = argparse.ArgumentParser(description='Return genome data for given accession IDs',
                                     epilog=example, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('accessionList', metavar='accessionID', nargs='+',
                        help='Accession IDs to return.')
    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')
    # The default will eventually be configurable with a config file.
    parser.add_argument('--cachedir', default=os.path.join(os.environ.get('HOME',''),'kSNP'),
                        help='Override the directory to use for caching data')
    # The default will eventually be configurable with a config file.
    parser.add_argument('--maxcachedays', type=int, default=90, help='Override the maximum age for cached data')
    parser.add_argument('--maxcachebytes', type=int, default=10000000000, help='Override the maximum size for the cache in bytes')
    parser.add_argument('--annotation', dest='responses', action='append_const', const=kcache.annotationData,
                        help='Output annotation data for the accession ID if available.')
    parser.add_argument('--genome', dest='responses', action='append_const', const=kcache.genomeData,
                        help='Output genome data for the accession ID if available.')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)


def getAccessionStatus(accessionList):
    # Relies on Entrez.email being set; currently being done in __main__, maybe
    # not the best place?

    # Variable to store Entrez records for these IDs
    records = {}
    # Flag to indicate we experienced an error retrieving data.
    searchError = False
    
    # esummary doesn't take a python list of accession IDs, unfortunately.
    accessionQuery = ", ".join(accessionList)

    logging.debug("Retrieving status for accession numbers: %s", accessionQuery)
    handle = Entrez.esummary(db="nuccore", id=accessionQuery, retmode="xml")
    parsedRecords = Entrez.parse(handle)
    try:
        for record in parsedRecords:
            records[record['AccessionVersion']] = record
            logging.info("Found entry for %s", record['AccessionVersion'])
            logging.debug("Entry data:\n%s", str(record))
                
    except StopIteration:
        # we got the last record.  All done.
        logging.info("No more accession information in status search.")
    except RuntimeError as errorMessage:
        if str(errorMessage).startswith("Invalid uid"):
            # In this case, we got NO data.  Need to go back around and try again
            # unless we only tried one ID and got this message.
            searchError = True
            records = {}
            logging.info("Unable to find entry with this ID.  Error: %s", errorMessage)
        else:
            raise

    if searchError and len(accessionList) > 1:
        # If we don't know what accessionList entry resulted in the error, try them
        # individually.
        # This seems a bit recursive... and it is, but the single entry check will
        # terminate the recursion after at most two levels.
        for accession in accessionList:
            records.update(getAccessionStatus([accession,]))

    #
    # Sanity checking
    #
    
    if len(accessionList) > len(records):
        # Since the accession number _is_ an identifier in nuccore, this
        # suggests something went wrong.  Not sure exception is right, but...
        logging.info("Not all accession IDs found in Entrez nuccore database!")

    for accession in accessionList:
        if accession not in records:
            logging.info("Asked for accession ID %s, but not in the response.  Possibly the request doesn't have a version number?", accession)

    return records

def retrieveAndCacheGenomeData(accessionList, cacheDir, maxCacheBytes, dataType=kcache.annotationData):
    # Note that this should be called by getGenomeDataThroughCache() because
    # that will ensure that the accession number used is the versioned one,
    # rather than the unversioned one.
    
    # This is where the annotation data will be stored (indexed by accession
    # number, containing an array of lines of text.
    resultData = {}
    
    for accession in accessionList:
        
        # Because efetch can return multiple records at once, it is actually
        # possible to pass a python list as id= value, and efetch will return
        # all of them.  However, too much parsing is required (any) to separate
        # and distinguish the records returned, so postponing that process.
        #
        # For the future, transition to xml return mode, use Entrez.parse, and
        # convert the resulting data structure into something much closer to
        # the structure we will use.  This should let us cache meaningfully
        # and reduce load on the server by making ony one request for multiple
        # annotations. - JN - TODO FIXME XXX
        handle = Entrez.efetch(db="nuccore", id=accession, rettype=dataType, retmode="text") # retmode?

        resultData[accession] = handle.read()

        handle.close()

        logging.info("Retrieved new %s data for %s", dataType, accession)

        cacheFile = None

        try:
            cacheFile = open(kcache.cacheFileName( accession, cacheDir, dataType),"w")
            cacheFile.write(resultData[accession])
            cacheFile.close()
            logging.info("Successfully cached %s into %s", accession, kcache.cacheFileName( accession, cacheDir, dataType))
            kcache.shrinkCache(cacheDir, maxCacheBytes)

        except OSError:
            if cacheFile is not None:
                cacheFile.close()
            logging.warning("WARNING: could not write cache file to cache directory for %s.  File permissions or a directory not existing could cause this.", accession)
            if os.access(kcache.cacheFileName(accession, cacheDir, dataType), os.W_OK):  # If we have access to the cache file, then
                # Remove the cache file because it is in an unknown state.
                os.remove(kcache.cacheFileName(accession, cacheDir))

    for accession in accessionList:
        if not accession in resultData:  # If the info for this accession number wasn't retrieved, then
            logging.info("WARNING: could not retrieve %s from NCBI; it is not available for annotation", accession)
            if os.access(kcache.cacheFileName(accession, cacheDir, dataType), os.W_OK):  # If we have access to the cache file, then
                # Remove the cache file because it is out of date.
                os.remove(kcache.cacheFileName(accession, cacheDir, dataType))

    return resultData


def findCachedGenomeData(accessionList, cacheDir, dataType=kcache.annotationData):
    # The actual annotation data retrieved from the cache
    resultData = {}

    # The status of the cache file, e.g. file changed date.
    resultStatus = {}

    cacheFile = None
    
    for accession in accessionList:
        try:
            cacheName = kcache.cacheFileName(accession, cacheDir, dataType=dataType)

            # Find the creation time of the cached file.  Later it may
            # make sense to see if there is relevant info in the cache
            # file itself.
            resultStatus[accession] = {
                'cacheTimeEpoch': os.stat(cacheName).st_ctime, }
            
            cacheFile = open(cacheName,"r")
            resultData[accession] = cacheFile.read()
            cacheFile.close()
            logging.debug("Found %s in the cache at %s", accession, cacheName)
            
        except OSError:
            # If there is an OS error reading the cached data, we delete
            # the cached item if we can and move on.  This would be the
            # case if the file didn't exist, or was not readable.
            if cacheFile is not None:
                cacheFile.close()
            if os.access(cacheName, os.W_OK):
                os.remove(cacheName)

            # Make sure we don't reference this accession number in our
            # output / returned data; we don't have it in our cache.
            resultData.pop(accession, None)
            resultStatus.pop(accession, None)
            logging.debug("Error retrieving cached value for %s", accession)

    return (resultStatus, resultData)


def cacheCurrent(annotation, databaseStatus, maxCacheDays):
    # Determine from the annotation data and the database status info
    # whether our cached data is current and usable, or whether we need to
    # retrieve updated original data.

    # Initial testing value: Everything functions, the cached data is never
    # returned - JN FIXME TODO XXX BUG

    # config values
    # Maximum amount of time before refreshing the cache anyway.
    maxCachedLifetime = datetime.timedelta(days=maxCacheDays)

    # Database status data
    lastChangedDateString = [ int(dateinfo) for dateinfo in databaseStatus['UpdateDate'].split('/') ]
    lastChangedDate = datetime.date(lastChangedDateString[0], lastChangedDateString[1], lastChangedDateString[2])
    canonicalName = databaseStatus['AccessionVersion']
    Replacement = databaseStatus['ReplacedBy']

    # Cached file data
    cachedDate = datetime.date.fromtimestamp(annotation['cacheTimeEpoch'])

    logging.debug("last changed: %s, canonical name: %s, replaced by: %s, cached file created date: %s", databaseStatus['UpdateDate'], canonicalName, Replacement, cachedDate)

    # Logic to determine if the cache is current
   
    return ( cachedDate > lastChangedDate and
             len(Replacement) == 0 and
             datetime.date.today() - cachedDate < maxCachedLifetime )

            
def getGenomeDataThroughCache(accessionList, cacheDir, maxCacheDays, maxCacheBytes, dataType=kcache.annotationData):
    # Passed a list of accession numbers, return annotation information for
    # each, making sure that 1) if the information must be retrieved from a
    # remote source, it is cached, and 2) if it has already been cached and
    # the cache is still valid, it is retrieved from the cache.

    # The list of accession numbers that need to be retrieved from source.
    accessionListToRetrieve = []

    databaseStatus = getAccessionStatus(accessionList)

    # We are getting the cached status of the official accession number which
    # includes the version.  Note this may be different from the list passed
    # to this function.
    (cachedStatus, cachedData) = findCachedGenomeData(databaseStatus.keys(), cacheDir, dataType=dataType)
    
    for accession in databaseStatus.keys():
        if accession in cachedData:
            if cacheCurrent(cachedStatus[accession],databaseStatus[accession], maxCacheDays):
                # If we have a cached annotation, and it is current, then skip
                # this retrieval.
                logging.info("Cached data for %s found and current", accession)
                continue
            
        # Otherwise we didn't have this cached, or it was out of date somehow.
        accessionListToRetrieve.append(accession)

    retrievedData = retrieveAndCacheGenomeData(accessionListToRetrieve, cacheDir, maxCacheBytes, dataType)

    for accession in accessionListToRetrieve:
        if not accession in retrievedData:
            # Remove the cached data (if any), because it is invalid.
            cachedData.pop(accession, None)
            logging.info("Previously cached data for %s, now invalid, but no current data for this accession ID", accession)


    # Return the cached annotations overwritten by any that were more recently
    # retrieved.
    for (key, value) in retrievedData.items():
        cachedData[key]=value
    return cachedData



####################################
####################################
###
###
### Main
###
###
####################################
####################################


if __name__ == "__main__":

    Entrez.email = "ksnp-dev@kissake.net"
    Entrez.tool = "kSNP4"
    outputData = {}

    options = parseCommandline()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s')

    logging.debug('Commandline options, as parsed: %s', str(options))

    if options.responses == None:
        responses = ['gb'] # Default to return annotations if not specified.
    else:
        responses = options.responses

    for dataType in responses:
        outputData[dataType] = getGenomeDataThroughCache(options.accessionList, options.cachedir, options.maxcachedays, options.maxcachebytes, dataType=dataType)
        
    for outputSection in outputData.keys():
        for (accession, output) in outputData[outputSection].items():
            logging.debug("Outputting %s for %s", outputSection, accession)
            sys.stdout.write(output)
