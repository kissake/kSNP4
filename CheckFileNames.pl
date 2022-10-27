#!/usr/bin/perl


use warnings;
use strict;

=begin
CheckFileNames cheks all of the file names in a kSNP# input file for violations of the name rules which are
File naming rules:
no spaces in name
no all-numeric names
no characters other than 0-9 A-Z a-z _ .
no more than one . in a file name

It then checks for duplicate genomeNames

Writes an output file listing the lines in the input file that viollate these rules

Usage: CheckFileNames infileName
=cut

my $infileName = $ARGV[0];
my @errors = ();
my $line = '';
my @temp = ();
my @temp2 = ();
my @temp3 = ();
my @temp4 = ();
my $fileName = '';
my $lineCount = 0;
my $errorString = '';
my $fn = ''; #part of $fileName before the extension
my $extension = ''; #part of file name after extension
my %genomes = (); #key = genomeID value = number of occurrences

#open the input file
open (INFILE, $infileName) or die "Can't open $infileName for reading. $!";

#process each line in the file
while ($line = <INFILE>) {
	chomp $line;
	$lineCount++;
	#print "$lineCount\n";
	@temp = split/\t/, $line;
	if (scalar @temp > 0) { # Ignore blank lines.
	    # Only if the line is not blank...
	    @temp2 = split/\//, $temp[0];	
	    $fileName = $temp2[$#temp2];

	    #test for more than one '.' in the file name
	    @temp3 = split /\./, $fileName;
	    #define the $fn
	    @temp4 = split /\./, $fileName;
	    $fn = $temp4[0];
	    if (scalar @temp3 > 2) {
		$errorString = "Line $lineCount:\t$fileName\n";
		push (@errors, $errorString);
	    }
	    elsif($fileName =~ /\:|\"|\'|\ |\?|\+|\)|\(|\*|\&|\^|\%|\$|\#|\@|\!|\;|\>|\</) { #check for spaces in names
		$errorString = "Line $lineCount:\t$fileName\n";
		push (@errors, $errorString);
	    }
	    #put each genomeID into %genomes 
	    if (exists $genomes{$temp[1]}) {
		$genomes{$temp[1]}++;
	    }
	    else {
		$genomes{$temp[1]} = 1;  
	    }
	}
}
close INFILE;

if (scalar keys (%genomes) > 0) {
	foreach my $genomeID(keys %genomes) {
		if ($genomes{$genomeID} == 1) {
			delete $genomes{$genomeID};
			}		
		}
	}


if (scalar @errors > 0 or  scalar keys(%genomes) > 0) {
	open (OUTFILE, ">NameErrors.txt") or die "Can't open NameErrors.txt for writing. $!";
	}

if (scalar @errors > 0) {	
	print "\n\nOne or more file names violates the rules for naming files.\n";
	print "The  file naming rules, followed by a list of the names that violate those rules\n";
	print "have been written to the file NameErrors.txt\n\nkSNP3 has been terminated.\n";
	
	print OUTFILE "One or more file names in the input file $infileName is illegal\n\n";
	print OUTFILE "The rules for naming files are:\n1.\tA file name may contian only one dot \(\'.\'\) character, that which separates the file ID from the extension.\n";
	print OUTFILE "\tEcoSME175.fasta is legal, EcoSME17.5.fasta is not\n\n";
	print OUTFILE "2.\tFile names may not contain any spaces.\n\tCF0123.fasta is legal, CF 0123.fasta is not\n\n";
	print OUTFILE "3.\tFile names may contain only the characters A-Z a-z 0-9 - . and _ \(the underscore character\).\n";
	print OUTFILE "\tForbidden characters include  : # + > < among others\n";
	print OUTFILE "\tEcoO157H7.fasta is legal, EcoO157:H7 is not.\n\n";
	print OUTFILE "Replace illegal characters with the underscore (_) if necessary for clarity.\n\n";
	print OUTFILE "The file names listed below, along with the line number in $infileName, are illegal\n";
	print OUTFILE "Be sure to correct illegal file names both on the file itself and in the input file $infileName.\n\n";
	print OUTFILE "After all file names have been corrected be sure to trash this file (NameErrors.txt).\n";
	print OUTFILE "If kSNP3 detects a file named NameErrors.txt it will terminate the run.\n\n";
	for(my $i = 0; $i < scalar @errors; $i++) {
		print OUTFILE $errors[$i];
		}
	}
	
if (scalar keys (%genomes) > 0) {
	print OUTFILE "\n\nThe following genome names occur multiple times.  Make each genomeID unique in the .in file.\n\n";
	foreach my $genomeID(keys %genomes) {
		if ($genomes{$genomeID} > 1) {
			print OUTFILE "$genomeID occurs $genomes{$genomeID} times.\n";
			}
		}
	print OUTFILE "\nBe sure to correct duplicate genome names in the input file $infileName.\n\n";
	print OUTFILE "After all file names have been corrected be sure to trash this file (NameErrors.txt).\n";
	print OUTFILE "If kSNP3 detects a file named NameErrors.txt it will terminate the run.\n\n";
	
	}
	
if (scalar @errors > 0 or scalar keys (%genomes) > 0) {
	close OUTFILE;
	}
