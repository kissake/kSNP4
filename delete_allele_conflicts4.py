#!/usr/bin/python3

import sys

def maskSNP(kmer,halfK):
    return kmer[:halfK] + "." + kmer[halfK+1:]


if __name__ == "__main__":

    # Assumptions:
    #  - Input is stdin, output is stdout (probably not true?)
    #  - Input is a list of kmers and their frequencies, but there is the off
    #    chance that some appear more than once (OOPS!)
    #  - Input need not be sorted.  If it is, no huge win, but if we _know_
    #    it is, we can change the algorithm a little to be more memory efficient.
    

    separator = " " # Default separator is a space character.
    newline = "\n" # Newline character for output.

    
    delete = {}  # Holds masks of kmers we are deleting.
    testing = {} # Holds masks and kmers for kmers we are considering keeping (no conflicts yet)
    
    line = sys.stdin.readline()

    # Determine k-mer length for this file from the first element.
    (kmer, freq) = line.split(" ",1)
    # We use this calculated value for clarity; it is (k-1) / 2.  We capture
    # the first halfK and last halfK characters from the k-mer when masking.
    halfK = int( (len(kmer) - 1) / 2 )

    while line:
        (kmer, freq) = line.split(" ",1)

        masked = maskSNP(kmer, halfK)
        if masked in delete.keys():
            pass # Nothing to do, we already know we're removing this one.
        elif masked in testing.keys():
            # This is probably a conflict, but we check (not sure input is unique):
            if kmer != testing[masked]:  # Otherwise, this is a dupe, and dupes are fine.
                # Now we need to remove this one
                delete[masked] = True
                del testing[masked] # Remove from kmers we are considering.
        else:
            testing[masked] = kmer

        line = sys.stdin.readline()

    # Found all of our conflicting k-mers, output the remaining ones:
    sys.stdout.write(newline.join(testing.values()))
                     
    #for kmer in testing.values():  # Note that we are extracting values instead of keys here.
    #    sys.stdout.write(kmer + newline)

