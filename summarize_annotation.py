#!/usr/bin/python2

import sys

# These hash tables are used as counters of unique items.  This works
# by setting the value of (e.g.) counters[x] to 1.  No matter how many
# times you do that, x shows up as a key only once, and y, which was
# never referenced, never shows up.  No sorting required.
counters={}
allLoci={}
headerOut=[]
statsOut=[]
tab="	"

def safeRatio( A, B ):
    if B == 0:
        return "inf"
    else:
        return A / B

# Omit header: 
headerline = sys.stdin.readline()
if not headerline.startswith("LocusNum"):
    print("Missing header!?!?")
    sys.exit(1)

for line in sys.stdin.readlines():
    print(line)
    ( locusId, context, nonSyn, rest ) = line.split(tab, 4) # Split on tabs
    if rest.find(tab):
        ( annotType, rest ) = rest.split(tab,2)
    else:
        annotType = rest

    allLoci[locusId] = 1 # One hash table entry for every locus found.
    annotationCounter = counters.setdefault(annotType, {})
    annotationCounter[locusId] = 1
    if nonSyn > 0:
        nonsynCounter = counters.setdefault("nonSyn", {})
        nonsynCounter[locusId] = 1 # Permits us to count non-synonymous SNPs.

unannotated = len(counters["UnannotatedRegion"])
nonsynonymous = len(counters["nonSyn"])
synonymous = len(allLoci) - unannotated - nonsynonymous 
annotated = len(allLoci) - unannotated

# Semicolons are rarely used in python.  Used here to make clear how the output is calculated.
headerOut.append("Num_NotAnnotatedRegion")	; statsOut.append( unannotated )
headerOut.append("AnnotatedNotProtein")		; statsOut.append( len(counters["NotProteinCoding"]) )
headerOut.append("Num_NonSynon")			; statsOut.append( nonsynonymous )
headerOut.append("Num_Synon")				; statsOut.append( synonymous )
headerOut.append("NS/S")					; statsOut.append( safeRatio(nonsynonymous, synonymous) )
headerOut.append("NSfractionOfAnnotated")	; statsOut.append( safeRatio(nonsynonymous, annotated) )
headerOut.append("NumLoci")					; statsOut.append( len(allLoci) )
headerOut.append("Num_InAnnotatedGenomes")	; statsOut.append( annotated )
headerOut.append("Num_NotInAnnotatedGenome"); statsOut.append( len(counters["NotInAnnotatedGenome"]) )

# Write output
print(tab.join(headerOut))
print(tab.join(statsOut))

