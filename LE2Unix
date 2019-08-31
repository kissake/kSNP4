#!/usr/bin/perl
=begin
Takes a list of paths and genome names as th input (same as kSNP3)
checks each file and returns the kind of line endings it has
Usage LE2Unix inputList
=cut
use warnings;
use strict;

my $listFile = $ARGV[0];
my $genomeID = ''; 
my $filePath = ''; #path to the file corresponding to the $genomeID
my $line = '';
my @temp = ();
my @Files = (); #col0 is fileID col1 is path
my $LE = '';
my $fileName = '';
my $pathToDirectory = '';
my @args = ();
my $macFlag = 'F';



#get the list of genomes and paths into @Files
open(INFILE, $listFile) or die "Can't open $listFile for reading. $!";
while ($line = <INFILE>) {
	chomp $line;
	@temp = split /\t/, $line;
	$genomeID = $temp[1];
	$filePath = $temp[0];
	push @Files, [$genomeID, $filePath];
	}
close INFILE;

#test each file to determine its file endings
for (my $i = 0; $i <scalar @Files; $i++) {# scalar @Files
	$LE = getLE($Files[$i][1]);
	@temp = split /\//, $Files[$i][1];
	$fileName = pop @temp;
	$pathToDirectory = join "\/", @temp;	
	if ($LE ne 'Unix'){
		print "The line endings of $fileName are $LE\n";
		}
	
	if ($LE eq 'Windows') {
		open(INFILE, $Files[$i][1]) or die "Can't open $Files[$i][1] for reading. $!";
		open (OUTFILE,'>temp.txt') or die "Can't open temp.txt for writing. $!";
		while ($line = <INFILE>) {
			$line =~ s/\r\n/\n/;
			print OUTFILE "$line";
			}
		close INFILE;
		close OUTFILE;
		rename ('temp.txt', $fileName);
		@args = ('mv', "$fileName", "$pathToDirectory");
		system(@args) == 0 or die "system @args failed: $!";
		}
		elsif($LE eq 'Mac') {
			print "\nWARNING!!! $fileName has Classic Mac line endings.\nYou should manually change the line endings to Unix!\n";	
			$macFlag = 'T';
			}
	}

if ($macFlag eq 'T') {
	print "\n\nThe files listed above have classic Mac line endings.\nYou should manually change those endings to Unix.\n\n";
	die;
	}


sub getLE
{
my $file = $_[0];
my $line = '';
my $lastChar = '';
my $secondLastChar = '';
my $LE = ''; #Windows, Unix or Mac

open (INFILE, $file) or die "Can't open $file for reading. $!";
$line = <INFILE>;
$lastChar = ord (substr($line, -1, 1));
$secondLastChar = ord(substr($line, -2, 1));
if($lastChar == 10 and $secondLastChar == 13) {
	$LE = 'Windows';
	}
	elsif ($lastChar == 10 and $secondLastChar != 13) {
		$LE = 'Unix';
		}
		elsif ($lastChar == 13) {
			$LE = 'Mac';
			}
close INFILE;
return $LE;
}