#!/usr/bin/python3

# Input:
# List of SNPs
# Tree information, including clades (either implicitly or explicitly)
# Tree file (newick format)
#
# Process the above information to identify the most likely root node
# for the phylogenic tree, and information about SNPs that occur between
# multiple phylogenic branches.
# 
# Output:
# Homoplasy groups
# ClusterInfo
# "NodeSigCounts" (file referenced on the commandline as last argument)
#

# Import standard Python libraries for argument parsing, interacting with the
# filesystem, and environment and time / date processing.
import argparse as argparse

# Permit use of stdout
import sys

# Help with debug data:
import logging as logging

# Phylogenic tree operations
import Bio.Phylo as Phylo


HomoplasticSNPsCountFile = 'COUNT_Homoplastic_SNPs'
HomoplasyGroupsFile = 'Homoplasy_groups'
ClusterInfoFile = 'ClusterInfo'




def parseCommandline(override=None):

    description = '''Generate phylogenetic trees based on SNP data and generated trees.  Determine most likely root
of the tree, and output informaiton about that tree and any identified  homoplasy groups.'''
    
    example = ''''''
    
    parser = argparse.ArgumentParser(description= description, epilog=example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('SNPsFile', metavar='SNPsFile', nargs=1,
                        help='The filename where discovered SNPs can be found.')

    parser.add_argument('NodeCladeStructures', metavar='NodeCladeStructures', nargs=1,
                        help='The filename where all possible clade structures for the nodes can be found.')

    parser.add_argument('Tree', metavar='Tree', nargs=1,
                        help='File containing the computed tree in Newick format.  The newly computed tree will use the same filename, with a .rerooted suffix.')

    parser.add_argument('OutputSigCounts', metavar='SigCounts', nargs=1,
                        help='Output filename for information about the tree nodes.')
    
    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')
    parser.add_argument('--info', action='store_true', help='Output status information on STDERR', default=True)
    parser.add_argument('--quiet', action='store_true', help='Silence debug and other messages. (warnings only)')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)


def inputClades(cladesFile):
    logging.debug("Starting to read nodes list")
    # Returns clades and IDsInNode
    logging.debug("Done reading nodes list")
    return None, None


def groupNodes(SNPsFile, clades, IDsInNodes):
    logging.debug("Starting to group nodes")
    # Returns core, group, gl, rootScore (was: num_map_to_node) and SNPsCount (# of SNPs)
    logging.debug("Done grouping nodes")
    return None, None, None, None, None

def findBestRoot(clades, rootScore):
    logging.debug("Starting to find best root")
    # Returns the best root, and that root's score (or the maximum number
    # of node loci associated with a given root)
    logging.debug("Done finding best root")
    return None, None

def writeRerootedTree(TreeFile, root):
    logging.debug("Starting to write rerooted tree")
    logging.warning("I don't know how to do the tree stuff...?!?")
    logging.debug("Done writing rerooted tree")
    exit(1)
    
def writeHomoplasticSNPCounts(HomoplasticSNPsCountFile, SNPsCount, maxRootScore):
    logging.debug("Starting to write homoplastic SNP counts")
    outputFile = open(HomoplasticSNPsCountFile, "w")
    # outputFile.write("Number_Homoplastic_SNPs: %s\n" % (SNPsCount - maxRootScore))
    outputFile.close()
    logging.debug("Done writing homoplastic SNP counts")

def writeNodeSigCounts(SigCountsFile, clades, root, group):
    logging.debug("Starting to write node sig counts")
    outputFile = open(SigCountsFile, "w")
    outputFile.close()
    logging.debug("Done writing node sig counts")
    return None
    
def writeHomoplasyGroups(HomoplasyGroupsFile, group, IDsInNode, root, clusters, nodeCount):
    logging.debug("Starting to write homoplastic SNP counts")
    outputFile = open(HomoplasyGroupsFile, "w")
    outputFile.close()
    logging.debug("Done writing homoplasy groups")


def writeClusterInfo(ClusterInfoFile, gl, clusters):
    logging.debug("Starting to write cluster info")
    outputFile = open(ClusterInfoFile, "w")
    outputFile.close()
    logging.debug("Done writing cluster info")


                     
if __name__ == "__main__":

    options = parseCommandline()

    if options.debug:
        logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%dT%H:%M:%S', level=logging.DEBUG)
    elif options.info:
        logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%dT%H:%M:%S', level=logging.INFO)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARN)
            

    logging.debug('Commandline options, as parsed: %s', str(options))

    # Input
    SNPsFile = options.SNPsFile[0]
    CladeStructuresFile = options.NodeCladeStructures[0]
    TreeFile = options.Tree[0]

    # Output
    SigCountsFile = options.OutputSigCounts[0]

          
    # Processing

    clades, IDsInNodes = inputClades(CladeStructuresFile)

    core, group, gl, rootScore, SNPsCount = groupNodes(SNPsFile, clades, IDsInNodes)
          
    root, maxScore = findBestRoot(clades, rootScore)
    
    writeHomoplasticSNPCounts(HomoplasticSNPsCountFile, SNPsCount, maxScore)

    clusters = writeNodeSigCounts(SigCountsFile, clades, root, group)

    nodeCount=0 # TODO FIXME it is actually:
    # nodeCount = len(clades[root]) # The number of distinct nodes within the tree.

    writeHomoplasyGroups(HomoplasyGroupsFile, group, IDsInNodes, root, clusters, nodeCount)

    writeClusterInfo(ClusterInfoFile, gl, clusters)

    writeRerootedTree(TreeFile, root)
    
