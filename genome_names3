#!/usr/bin/perl
#v3.1
# 
# Example: genome_names in_list annotated_genomes
#
# Prints to the file named in the second argument the genome names that are in you input list of genomes. 
#
use warnings;
use strict;

my $infile     = $ARGV[0];
my $outfile = $ARGV[1];
my $line = '';
my @temp = ();

open (INFILE, $infile) or die "Can't open $infile for reading. $!";
open (OUTFILE, '>'.$outfile) or die "Can't open $outfile for writing. $!";

while ($line = <INFILE>) {
	chomp $line;
	@temp = split /\t/, $line;
	print OUTFILE "$temp[1]\n";
	}
close INFILE;
close OUTFILE;