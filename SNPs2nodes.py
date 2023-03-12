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

# Regular expressions for input parsing ... Not a huge deal as long
# as our files aren't insane.

HomoplasticSNPsCountFile = 'COUNT_Homoplastic_SNPs'
HomoplasyGroupsFile = 'Homoplasy_groups'
ClusterInfoFile = 'ClusterInfo'

'''
Template for buffered reading of files:
    # Setup for buffered reading of input file
    bufferSize = 32 * 1024 # Make this a global, or commandline argument?
    inputFile = open(cladesFile, "r")
    nextLines = inputFile.readlines(bufferSize)

    # Loop setup.
    
    while nextLines != []:
        # This line should be all alone ; no other lines indented the same
        for line in nextLines:
        # This line should be all alone ; no other lines indented the same

            # Inner per-line loop.

        # This line should be all alone ; no other lines indented the same
        nextLines = inputFile.readlines(bufferSize)
        # This line should be all alone ; no other lines indented the same

    inputFile.close()

'''




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

    # Setup for buffered reading of input file
    bufferSize = 32 * 1024 # Make this a global, or commandline argument?
    inputFile = open(cladesFile, "r")
    nextLines = inputFile.readlines(bufferSize)

    # A dict to translate from genome names to corresponding bit values.
    genomeBits = {}
    maxGenomeBit = 1 # Very first genome found is 1, then increase by factor of 2

    clades = {}
    IDsInNode = {}
    
    currRoot = None
    currNode = None
    genomeList = 0
    
    while nextLines != []:
        for line in nextLines:
            # Valid lines either have:
            #   4 elements ("root:", root node, "node:", sub node),
            #   one element (leaf node), or
            #   No elements (empty line, end of clade)
            elements = line.strip().split()
            # logging.debug("Elements: %s", elements)
            if len(elements) == 0:
                # End of clade, looking for a new one.
                clades.setdefault(currRoot, {})[currNode] = genomeList
                IDsInNode.setdefault(genomeList, {})[currRoot] = currNode
                # logging.debug("End of clade w/ root: %s and node: %s.  Genomes: %s", currRoot, currNode, genomeList)
                # The below is immediately overwritten by start-of-clade behavior.
                currRoot = None
                currNode = None
                genomeList = 0
                
            elif len(elements) == 4:
                # Start of clade:
                if elements[0] == "root:" and elements[2] == "node:":
                    currRoot = elements[1]
                    currNode = elements[3]
                    genomeList = 0 # Empty bit field.
                    # logging.debug("Start of clade w/ root: %s and node: %s.", currRoot, currNode)
                                    
                else:
                    logging.warning("Syntax error reading clades from %s.  Exiting.", cladesFile)
                    exit(1)
            else:
                # Identify the bit corresponding to the current genome (or add this genome to the list
                # if it isn't there yet)
                thisGenome = elements[0]
                thisBit = genomeBits.setdefault(thisGenome, maxGenomeBit)
                if thisBit == maxGenomeBit:
                    # logging.debug("Geome: %s => %s", thisGenome, thisBit)
                    maxGenomeBit = maxGenomeBit * 2
            
                # Set the bit for this genome in this clade.
                genomeList = genomeList | thisBit

        nextLines = inputFile.readlines(bufferSize)

    # Handle the very last clade we reviewed, in case there was no terminal bare newline in the file:
    # End of clade
    clades.setdefault(currRoot, {})[currNode] = genomeList
    IDsInNode.setdefault(genomeList, {})[currRoot] = currNode
     
    inputFile.close()

    logging.debug("Length(clades): %s", len(clades.keys()))
    logging.debug("Length(IDsInNode): %s", len(IDsInNode.keys()))
    
    logging.debug("Done reading nodes list")
    return clades, IDsInNode, genomeBits


def groupNodes(SNPsFile, clades, idsInNodes, genomeBits):
    logging.debug("Starting to group nodes")

    ### Output variables:
    core = {}
    group = {}
    gl = {}
    rootScore = {}
    # An overall count (probably equal to greatest SNP loci ID)
    lociAddedCount = 0

    ### Symbols
    
    # Loop setup.
    tab = '\t'

    # Column identifiers
    lociID = 0 # The numeric ID used for this loci.
    center = 2
    loci = 1
    position = 3 # Maybe we don't need this?
    genome = 4

    # This should set every genome bit.  If python didn't have unlimited
    # width integers, this would probably overflow quickly.
    allGenomes = sum(genomeBits.values())


    ### Internal state
    
    # Ensure these variables have scope to maintain state by declaring
    # them outside of the loop.
    currentLocus = {}
    lastSnpLocus = None
    allVariants = {}
 
    # Setup for buffered reading of input file
    bufferSize = 32 * 1024 # Make this a global, or commandline argument?
    inputFile = open(SNPsFile, "r")
    nextLines = inputFile.readlines(bufferSize)

    while nextLines != []:
        for SNP in nextLines:

            # logging.debug("SNP: %s",SNP)
            # Inner per-line loop.
            snpData = SNP.strip().split(tab)
            if len(snpData) == 1:
                # Empty line
                continue
            snpCenter = snpData[center]
            snpLocus = snpData[loci]
            snpLocusID = snpData[lociID]
            snpPosition = snpData[position]
            snpGenome = snpData[genome]
            thisBit = genomeBits[snpGenome]
            
                
            if snpLocus != lastSnpLocus:
                if lastSnpLocus is None:
                    # This is actually our first one.
                    lastSnpLocus = snpLocus
                    continue
                if ( ( lociAddedCount % 10000 ) == 0 ):
                    logging.debug("Processed %s SNPs, currently on %s", lociAddedCount, lastSnpLocus)
                    
                ##### Execute for every loci in SNPs_all:
                # This code needs to be repeated at the very end of the loop
                # to ensure we hit the last SNP.
                lociAddedCount = lociAddedCount + 1

                # If this SNP appears in all genomes, it is a core genome; otherwise
                # it is False (not) (assigning the result of a boolean operation isn't
                # super-common, but is legit.)
                core[lastSnpLocus] = sum(allVariants.values()) == allGenomes

                # Create group and gl maps between genomes and locus ID:
                matchingRoots = {}
                for genomeSet in allVariants.values():
                    group.setdefault(genomeSet,{})[lastSnpLocus] = True
                    gl.setdefault(lastSnpLocus,{})[genomeSet] = True

                    # Identify roots that have at least one clade that
                    # matches at least one variant of this locus.
                    for root in idsInNodes.get(genomeSet,{}).keys():
                        matchingRoots[root] = True

                for root in matchingRoots.keys():
                    rootScore[root] = rootScore.get(root,0) + 1
                    
                
                ##### Reset for the current (new) locus.
                currentLocus = {}
                allVariants = {
                    'A': 0,
                    'C': 0,
                    'G': 0,
                    'T': 0,
                    '-': 0 # Not sure we keep this one...?
                }
                lastSnpLocus = snpLocus

            ##### Handle the current (new) locus
            currentLocus.setdefault(thisBit, {})[snpPosition] = snpCenter
            if snpCenter == '':
                snpCenter = '-'
            allVariants[snpCenter] = allVariants.get(snpCenter,0) | thisBit

                
        nextLines = inputFile.readlines(bufferSize)


    ##### Execute for every loci in SNPs_all ; processing the very last one
    # This code needs to be repeated at the very end of the loop
    # to ensure we hit the last SNP.
    lociAddedCount = lociAddedCount + 1

    # If this SNP appears in all genomes, it is a core genome; otherwise
    # it is False (not) (assigning the result of a boolean operation isn't
    # super-common, but is legit.)
    core[lastSnpLocus] = sum(allVariants.values()) == allGenomes

    # Create group and gl maps between genomes and locus ID:
    matchingRoots = {}
    for genomeSet in allVariants.values():
        group.setdefault(genomeSet,{})[lastSnpLocus] = True
        gl.setdefault(lastSnpLocus,{})[genomeSet] = True

        # Identify roots that have at least one clade that
        # matches at least one variant of this locus.
        for root in idsInNodes.get(genomeSet,{}).keys():
            matchingRoots[root] = True

    for root in matchingRoots.keys():
        rootScore[root] = rootScore.get(root,0) + 1

    ##### And now we are done processing the input file.


    inputFile.close()

    logging.debug("Done grouping nodes")

    logging.debug("Core SNPs: %s", len([ value for value in core.values() if value ]))

    logging.debug("Groups: %s", len(group.keys()))

    logging.debug("Loci pointing to groups: %s", len(gl.keys()))

    logging.debug("Roots: %s", len(rootScore.keys()))

    logging.debug("Loci: %s", lociAddedCount)
                  
    
    # Returns core, group, gl, rootScore (was: num_map_to_node) and SNPsCount (# of SNPs)
    return core, group, gl, rootScore, lociAddedCount

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

    clades, IDsInNodes, genomeBits = inputClades(CladeStructuresFile)

    core, group, gl, rootScore, SNPsCount = groupNodes(SNPsFile, clades, IDsInNodes, genomeBits)
          
    root, maxScore = findBestRoot(clades, rootScore)
    
    writeHomoplasticSNPCounts(HomoplasticSNPsCountFile, SNPsCount, maxScore)

    clusters = writeNodeSigCounts(SigCountsFile, clades, root, group)

    nodeCount=0 # TODO FIXME it is actually:
    # nodeCount = len(clades[root]) # The number of distinct nodes within the tree.

    writeHomoplasyGroups(HomoplasyGroupsFile, group, IDsInNodes, root, clusters, nodeCount)

    writeClusterInfo(ClusterInfoFile, gl, clusters)

    writeRerootedTree(TreeFile, root)
    
