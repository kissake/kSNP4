#!/bin/bash -evx

# For debugging, use below as first line.
#!/bin/bash -vx

# Examples are default in Example/Example<#>/...

# Since this script will be in the source, default to making sure that the
# binary dependencies (jellyfish, mummer, consense, etc.) are available:

echo "Checking dependencies:"
[ -x jellyfish ] 	&& echo "... Jellyfish found in current working directory (cwd)" 	|| echo "Jellyfish missing"
[ -x mummer ] 		&& echo "... mummer found in cwd" 					|| echo "mummer missing"
[ -x consense ] 	&& echo "... consense found in cwd" 				|| echo "consense missing"
[ -x FastTreeMP ] 	&& echo "... FastTreeMP found in cwd" 				|| echo "FastTreeMP missing"
[ -x parsimonator ] 	&& echo "... parsimonator found in cwd" 				|| echo "parsimonator missing"
[ -x kSNP4 ] 		&& echo "... kSNP4 found in cwd" 					|| echo "kSNP4 missing"

MISSINGEXAMPLES=0
echo "Checking for required example files:"
[ -d Examples/ ] || MISSINGEXAMPLES=1 
[ -d Examples/Example2 ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/annotated_genomes ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/in_list2 ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/in_list ] || MISSINGEXAMPLES=1 
[ -d Examples/Example2/Genomes ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/Genomes/Vc63-93_MO45.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/Genomes/VcO1-ElTorN16961.fa ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/Genomes/VcMS6.fa ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/Genomes/Vc523-80.fa ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/Genomes/VcLMA3984-4.fa ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/Genomes/ERR579925.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example2/CommandLines.txt ] || MISSINGEXAMPLES=1 
[ -d Examples/Example1 ] || MISSINGEXAMPLES=1 
[ -d Examples/Example1/Genomes ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_PE-0_0155.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_NJ-60.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_BeAr436087.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_FL93-939.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_ref_gi21218484.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_Georgia97.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_PE6.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_Florida91-4697.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_PE-3_0815.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/Genomes/EEE_NorthAmerican_antigenic_variety.fasta ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/ExampleRuns/in_list ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/ExampleRuns/annotated_genomes ] || MISSINGEXAMPLES=1 
[ -f Examples/Example1/CommandLines.txt ] || MISSINGEXAMPLES=1 

if [ ${MISSINGEXAMPLES} -gt 0 ] 
then
	echo "Missing critical files under Examples/; are you sure you have unpacked the Examples?" 
	exit 1
else
	echo "... Found examples."
fi


# Typically the in_list file must be modified to run the examples.

# The output directory and log output locations can be specified.  There may 
# be hidden path dependencies such that the output directory must be a sub-
# directory of the current.

# Current working directory:
CURRENT=`pwd`

# Directory for caching
CACHE="${HOME}/kSNP/"

# Create our test directory (temporary directory)
TEMPDIR=`mktemp -d -p "${CURRENT}"`
INLIST1="${TEMPDIR}/in_list1"
ANNDGENOMES1="${CURRENT}/Examples/Example1/ExampleRuns/annotated_genomes"

RUN11DIR="${TEMPDIR}/Run11"
RUN11LOG="${TEMPDIR}/Run11Log"
RUN11ERRORLOG="${TEMPDIR}/Run11ErrorLog"

RUN12DIR="${TEMPDIR}/Run12"
RUN12LOG="${TEMPDIR}/Run12Log"
RUN12ERRORLOG="${TEMPDIR}/Run12ErrorLog"

RUN13DIR="${TEMPDIR}/Run13"
RUN13LOG="${TEMPDIR}/Run13Log"
RUN13ERRORLOG="${TEMPDIR}/Run13ErrorLog"

INLIST21="${TEMPDIR}/in_list21"
INLIST22="${TEMPDIR}/in_list22"
ANNDGENOMES2="${CURRENT}/Examples/Example2/annotated_genomes"

RUN21DIR="${TEMPDIR}/Run21"
RUN21LOG="${TEMPDIR}/Run21Log"
RUN21ERRORLOG="${TEMPDIR}/Run21ErrorLog"

RUN22DIR="${TEMPDIR}/Run22"
RUN22LOG="${TEMPDIR}/Run22Log"
RUN22ERRORLOG="${TEMPDIR}/Run22ErrorLog"

# Copy first in_list into place
sed -e "s/^.*\/Genomes/Genomes/"  < Examples/Example1/ExampleRuns/in_list | while read LINE ; do echo "${CURRENT}/Examples/Example1/${LINE}" ; done > "${INLIST1}"
echo "Generated in_list for Example1:"
cat "${INLIST1}"

### # Example 1
### Run1:
### kSNP4 -in in_list -outdir Run1 -k 13 | tee Run1Log
### 
### Run2:
### kSNP4 -in in_list -outdir Run2 -k 13 -annotate annotated_genomes | tee Run2Log
### 
### Run3:
### kSNP4 -in in_list -outdir Run3 -k 13 -annotate annotated_genomes -ML -NJ -vcf -core -min_frac 0.75 | tee Run3Log
### 

kSNP4 -cachedir "${CACHE}" -debug -in "${INLIST1}" -outdir "${RUN11DIR}" -k 13 2>&1 | tee "${RUN11LOG}" 

kSNP4 -cachedir "${CACHE}" -debug -in "${INLIST1}" -outdir "${RUN12DIR}" -k 13 -annotate "${ANNDGENOMES1}" 2>&1 | tee "${RUN12LOG}" 

kSNP4 -cachedir "${CACHE}" -debug -in "${INLIST1}" -outdir "${RUN13DIR}" -k 13 -annotate "${ANNDGENOMES1}" -ML -NJ -vcf -core -min_frac 0.75 2>&1 | tee "${RUN13LOG}" 

### # Example2
### 
### Run1:
### kSNP4 -in in_list -outdir Run1 -k 19 -annotate annotated_genomes | tee Run1Log
### 
### Run2: Note: replace ### with the path to SNPs_all in the Run1 directory!  If any directory names in the path contain spaces be sure to enclose those directory names in quotes (" ").
### 
### kSNP4 -in in_list2 -outdir Run2 -k 19 -annotate annotated_genomes -SNPs_all ###Examples/Example2/ExampleRuns/Run1/SNPs_all### | tee Run2Log

# Copy in_list into place
sed -e "s/^.*\/Genomes/Genomes/"  < Examples/Example2/in_list | while read LINE ; do echo "${CURRENT}/Examples/Example2/${LINE}" ; done > "${INLIST21}"
echo "Generated in_list for example2, run 1:"
cat "${INLIST21}"

# Copy in_list into place
sed -e "s/^.*\/Genomes/Genomes/"  < Examples/Example2/in_list2 | while read LINE ; do echo "${CURRENT}/Examples/Example2/${LINE}" ; done > "${INLIST22}"
echo "Generated in_list for example2, run 2:"
cat "${INLIST22}"


{ # The second has a dependency on the first. 
	kSNP4 -cachedir "${CACHE}" -debug -in "${INLIST21}" -outdir "${RUN21DIR}" -k 19 -annotate "${ANNDGENOMES2}" 2>&1 | tee "${RUN21LOG}"

	kSNP4 -cachedir "${CACHE}" -debug -in "${INLIST22}" -outdir "${RUN22DIR}" -k 19 -annotate "${ANNDGENOMES2}" -SNPs_all "${RUN21DIR}/SNPs_all" 2>&1 | tee "${RUN22LOG}" 
} 

echo "Waiting for all examples to complete."
wait

echo "Exiting for review: ${TEMPDIR}"
exit

### This code never reached.  Temporary directory must be cleaned up manually.

# Clean up temporary test directory
rm -r "${TEMPDIR}"

