#!/usr/bin/python3

# Import the Entrez database API from BioPython
import Bio.Entrez as Entrez

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


# The override option to this function can be used to cause the parsing of a provided list
# instead of the programs passed options, for testing.
def parseCommandline(override=None):
    parser = argparse.ArgumentParser(description='Return genome accession data for given accession IDs')
    parser.add_argument('accessionList', metavar='accessionID', nargs='+',
                        help='Accession IDs to return.')
    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')
    # The default will eventually be configurable with a config file.
    parser.add_argument('--cachedir', default=os.environ.get('HOME','') + '/kSNP/GbkCache/', help='Override the directory to use for caching data')
    # The default will eventually be configurable with a config file.
    parser.add_argument('--maxcachedays', type=int, default=90, help='Override the maximum age for cached data')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)
    
    
def cacheFileName(accession, cacheDir):
    # This shouldn't be just a string concatenation; what about directory separators?
    return cacheDir + str(accession)

def getAccessionStatus(accessionList):
    # Relies on Entrez.email being set; currently being done in __main__, maybe
    # not the best place.

    # Variable to store Entrez records for these IDs
    records = {}
    
    # esummary doesn't take a python list of accession IDs, unfortunately.
    accessionQuery = ", ".join(accessionList)

    handle = Entrez.esummary(db="nuccore", id=accessionQuery, retmode="xml")
    for record in Entrez.parse(handle):
        # Add the record to our dict, indexed by accession number + version.
        records[record['AccessionVersion']] = record

    if len(accessionList) > len(records):
        # Since the accession number _is_ an identifier in nuccore, this
        # suggests something went wrong.  Not sure exception is right, but...
        logging.warning("WARNING: not all accession IDs found in Entrez nuccore database!")
        raise KeyError

    for accession in accessionList:
        if accession not in records:
            logging.info("WARNING: Asked for accession ID %s, but not in the response.  Possibly the request doesn't have a version number?", accession)

    return records

def retrieveAndCacheAnnotations(accessionList, cacheDir):
    # Note that this should be called by getAnnotationsThroughCache() because
    # that will ensure that the accession number used is the versioned one,
    # rather than the unversioned one.
    
    # This is where the annotation data will be stored (indexed by accession
    # number, containing an array of lines of text.
    annotationData = {}
    
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
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") # retmode?

        annotationData[accession] = handle.read()

        handle.close()

        logging.info("Retrieved new data for %s", accession)

        cacheFile = None

        try:
            cacheFile = open(cacheFileName( accession, cacheDir),"w")
            cacheFile.write(annotationData[accession])
            cacheFile.close()
            logging.info("Successfully cached %s", accession)

        except OSError:
            if cacheFile is not None:
                cacheFile.close()
            logging.warning("WARNING: could not write cache file to cache directory for %s.  File permissions or a directory not existing could cause this.", accession)
            if os.access(cacheFileName(accession, cacheDir), os.W_OK):  # If we have access to the cache file, then
                # Remove the cache file because it is in an unknown state.
                os.remove(cacheFileName(accession, cacheDir))
            
    '''
    for record in Entrez.parse(handle):
        annotationData[record.name] = record
        cacheFiles[record.name].write(record)

    for file in cacheFiles.values():
        file.close()
    '''

    for accession in accessionList:
        if not accession in annotationData:  # If the info for this accession number wasn't retrieved, then
            logging.info("WARNING: could not retrieve %s from NCBI; it is not available for annotation", accession)
            if os.access(cacheFileName(accession, cacheDir), os.W_OK):  # If we have access to the cache file, then
                # Remove the cache file because it is out of date.
                os.remove(cacheFileName(accession, cacheDir))

    return annotationData


def findCachedAnnotation(accessionList, cacheDir):
    # The actual annotation data retrieved from the cache
    annotationData = {}

    # The status of the cache file, e.g. file changed date.
    annotationStatus = {}

    cacheFile = None
    
    for accession in accessionList:
        try:
            cacheName = cacheFileName(accession, cacheDir)

            # Find the creation time of the cached file.  Later it may
            # make sense to see if there is relevant info in the cache
            # file itself.
            annotationStatus[accession] = {
                'cacheTimeEpoch': os.stat(cacheName).st_ctime, }
            
            cacheFile = open(cacheName,"r")
            annotationData[accession] = cacheFile.read()
            cacheFile.close()
            logging.debug("Found %s in the cache", accession)
            
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
            annotationData.pop(accession, None)
            annotationStatus.pop(accession, None)
            logging.debug("Error retrieving cached value for %s", accession)

    return (annotationStatus, annotationData)


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

            
def getAnnotationsThroughCache(accessionList, cacheDir, maxCacheDays):
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
    (cachedStatus, cachedAnnotations) = findCachedAnnotation(databaseStatus.keys(), cacheDir)
    
    for accession in databaseStatus.keys():
        if accession in cachedAnnotations:
            if cacheCurrent(cachedStatus[accession],databaseStatus[accession], maxCacheDays):
                # If we have a cached annotation, and it is current, then skip
                # this retrieval.
                logging.info("Cached data for %s found and current", accession)
                continue
            
        # Otherwise we didn't have this cached, or it was out of date somehow.
        accessionListToRetrieve.append(accession)

    retrievedAnnotations = retrieveAndCacheAnnotations(accessionListToRetrieve, cacheDir)

    for accession in accessionListToRetrieve:
        if not accession in retrievedAnnotations:
            # Remove the cached data (if any), because it is invalid.
            cachedAnnotations.pop(accession, None)
            logging.info("Previously cached data for %s, now invalid, but no current data for this accession ID", accession)


    # Return the cached annotations overwritten by any that were more recently
    # retrieved.
    return cachedAnnotations | retrievedAnnotations
    

if __name__ == "__main__":

    Entrez.email="ksnp-dev@kissake.net"

    options = parseCommandline()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s')
    
    annotations = getAnnotationsThroughCache(options.accessionList, options.cachedir, options.maxcachedays)

    for annotation in annotations.values():
        sys.stdout.write(annotation)
    
