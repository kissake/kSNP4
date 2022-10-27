#!/usr/bin/perl -w
#v3.0

#rename_from_table.pl file_for_renaming fasta.table file.renamed
no warnings 'deprecated';

my $infile=$ARGV[0];
my $fasta_table=$ARGV[1];
my $out=$ARGV[2];
if (!defined $out ) {
    $out=$infile.".renamed";
}

open IN,"$infile";
@lines=<IN>;
close IN;
open IN,"$fasta_table";
@table=<IN>;
close IN;
open OUT,">$out";

my %names=();
my $prefix;
foreach my $entry (@table) {
    $entry =~ s/PREPROCESSED//;
    if ($entry =~ /(\w+)\s+(?:\>|)(.*)/) {
	$idx=$1;
	$name=$2;
	$names{$idx}=$name;
	if ($idx=~/(.*?)[0-9]+/) {
	    $prefix=$1;
	    #print "$prefix\n";
	}
    }
}
foreach $line (@lines) {
    #print "initial:  $line";
    if ($line =~ /($prefix[0-9]+)/) { 
	$idx=$1;
#    foreach my $idx (keys %names) {
	$line =~ s/$idx(\:|\s)/$names{$idx}$1/;
    }
    #$line =~ s/x\:/\:/g;
    print OUT "$line";
}
close OUT;
