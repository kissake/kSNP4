#!/usr/bin/python3

import Bio as Bio
import configparser as configparser
import argparse as argparse
import sys as sys
import os as os

configFileName = '.kSNPconfig'

configItems = {
    'CACHEDIR': {
        'help': 'This is the directory that kSNP will use to cache any time- or network-intensive work it performs.  Enter "None" to disable caching.',
        'prompt': 'Cache directory',
        'default': os.path.join(os.environ.get('HOME', None), 'kSNP'),
    },
    'CACHEBYTES': {
        'help': 'To avoid consuming all of the available disk space, kSNP will limit its use of the cache directory to the specified number of bytes',
        'prompt': 'Maximum cache size (bytes)',
        'default': '10000000000',
    },
    'MAXCACHEDAYS': {
        'help': 'For data that cannot be validated as consistent with the source, but which we want to cache if possible, how long should we retain the data before automatically freshening it from the source?',
        'prompt': 'Maximum cache staleness (days)',
        'default': '90'
    },
    'NUMCPUS': {
        'help': 'To make most efficient use of resources, kSNP attempts to use all of the available CPUs by default.  You can override that behavior by specifying the preferred number of processes to run concurrently here.',
        'prompt': 'Number of CPUs to use',
        'default': '0', # this means automatic
    },
    'VERBOSITY': {
        'help': 'The default amount of output you want to see from kSNP.  This can be used to make kSNP completely silent except if there is an error (Verbosity = 0), or exceedingly noisy, sharing details of every command issued (Verbosity = 2).',
        'prompt': 'Verbosity level (0, 1, 2)',
        'default': '1', # This means progress meter, but not debug.
    },
}

# For future consideration, the following possible config items:
# - Genome search directory (to avoid having to use full paths in input files)
# - Log file creation (so user doesn't have to redirect output)
# - Specific paths for genome or annotation data or other caches.
# - Whether to compress cached data.
# - Whether to cache SNPs and other output data (database)



def parseCmdline(alternateArgs=None):
    pass

def parseConfig(configfile=None):
    pass

def promptWithDefault(prompt, default=None):
    # Prompts the user for a value, offering a default if appropriate.
    if default is not None:
        response = input("%s [%s]: " % ( prompt, default ))
    else:
        response = input("%s: " % ( prompt ))

    if len(response) == 0:
        return default
    else:
        return response



def promptForUpdateConfig(oldConfiguration=None):
    newConfig = configparser.ConfigParser()
    newConfig['main'] = {}

    for configItem in configItems.keys():
        if oldConfiguration is None:
            newConfig['main'][configItem] = promptWithDefault(configItems[configItem]['prompt'], configItems[configItem]['default'])

        else:
            newConfig['main'][configItem] = promptWithDefault(configItems[configItem]['prompt'], oldConfiguration['main'][configItem])

    if newConfig['main']['CACHEDIR'] == 'None':
        # Conflict between there being a default and a potentially desirable setting being the empty string.
        newConfig['main']['CACHEDIR'] = ''

    return newConfig
        
        
if __name__ == '__main__':
    configFile = os.path.join(os.environ.get('HOME'),configFileName)
    if os.access(configFile,  os.F_OK):
        config = configparser.ConfigParser()
        config.read(configFile)
    else:
        config = None
    
    updatedConfig = promptForUpdateConfig(config)
    writeConfig = open(configFile,'w')
    updatedConfig.write(writeConfig)
    writeConfig.close()
