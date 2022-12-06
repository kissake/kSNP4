#!/bin/bash

### INPUTS
# Matrix file
# Job name
# Number of trees
# Number of concurrent jobs

### OUTPUTS
# "intree" (input to consense)

### Internal:
# Need a source of seed info, ideally such that the same seeds are generated
# for the same number of trees and number of jobs.
# Answer:  Bash's RANDOM variable.

# Seed bash random number generator
RANDOM=1234
PARSIMONATOR=/home/jnisbet/Documents/Development/ksnp/parallel_trees/kSNP4/binaries/parsimonator

MATRIX="${1}"
JOBID="${2}"
TREES="${3}"
CPUS="${4}"
OUTPUT="${5}"

let TREESPERCPU="${TREES}"/"${CPUS}"

TREESTODO="${TREES}"
CPUSUSED=0

while [ "${CPUSUSED}" -lt "${CPUS}" ]
do
    let CPUSUSED="${CPUSUSED}"+1
    let TREESTODO="${TREESTODO}"-"${TREESPERCPU}"
    if [ "${TREESTODO}" -lt "${TREESPERCPU}" ]
    then
	# If the number of trees is not evenly divisible
	# by the number of CPUs, our last parsimonator run
	# won't need to make as many trees.
	TREESPERCPU="${TREESTODO}"
    fi
    "${PARSIMONATOR}" -s "${MATRIX}" -n "${JOBID}.${CPUSUSED}" -N "${TREESPERCPU}" -p "${RANDOM}" 2>&1 > /dev/null &
done

# Wait for all of the parsimonator jobs to finish
wait

# Output all parsimony trees corresponding to the best (shortest) length to output file
TEMPINFO=treescores

grep "^Parsimony tree" RAxML_info.${JOBID}* | sort -k6 -n | awk '{print $6 " " $14}' > "${TEMPINFO}"

TOPSCORE=`head -1 "${TEMPINFO}" | cut -d\  -f 1`

cat "${TEMPINFO}" | while read SCORE NAME
do
    if [ "${SCORE}" -eq "${TOPSCORE}" ]
    then
	cat "${NAME}" >> "${OUTPUT}"
    else
	break
    fi
done
