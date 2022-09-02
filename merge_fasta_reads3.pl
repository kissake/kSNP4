#!/usr/bin/perl
#v3.1

#merge_fasta_contigs.pl $fasta > $fasta.merged
#Program to take a fasta file and merges all the entries into one, each entry separated by 1 N and strings of N's replaced by a single N. 

# Example: merge_fasta_contigs.pl contig_fasta_input > merged_fasta_output
# fasta header will be name of the input file plus " merged"
no warnings 'deprecated';

my $fasta = $ARGV[0];
#print "fasta=$fasta\n";

#my $fileName_noPath=`basename $fasta`;
if ($fasta =~ /.*\/(.*?)$/) {
    $fileName_noPath=$1;
}

chomp $fileName_noPath;
$fileName_noPath =~ s/\s+//;

#print "$fileName_noPath\n";
my $filler="N";
open IN,"$fasta" or die "Can't open $fasta for reading: $!\n";
my $id="9999999nonsense";


while ($line = <IN>) {
    chomp $line;  
    if ($line =~ /^>(.*)/ && $id ne "9999999nonsense") {
	print "$filler";
    }
    if ($line =~ /^>(.*)/ && $id eq "9999999nonsense") {
	#$id = $1." merged";
	$id = $fileName_noPath."  merged";
	print ">$id\n";
    }

    if ($line !~ /^>(.*)/) {
	$line =~ s/N+/N/g;
	print "$line";
    }
 
}   
print "\n";