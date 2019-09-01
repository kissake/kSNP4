#!/usr/bin/perl -w
#v3.0

#rename_from_table.pl file_for_renaming fasta.table file.renamed
no warnings 'deprecated';

# all_SNPs_sorted_labelLoci fileName2genomeName SNPs_all

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
    # NOT 100% clear what this regex does.  TODO review and re-evaluate - JN
    if ($entry =~ /(\w+)\s+(?:\>|)(.*)/) {
	# Apparent thrust is to turn fileName2genomeName into a hash table to
	# permit easy translation in the second loop. - JN
	$idx=$1;
	$name=$2;
	$names{$idx}=$name;
	# TODO: This is run every loop; can probably be just run for the first
	# entry?  No performance impact, as these files are v. small. - JN
	if ($idx=~/(.*?)[0-9]+/) {
	    $prefix=$1;
	    #print "$prefix\n";
	}
    }
}
foreach $line (@lines) {
    #print "initial:  $line";
    # IF and only if the line has a filename prefix from above, replace filename
    # with genome name from file processed above, with the filename appended(?) - JN
    if ($line =~ /($prefix[0-9]+)/) { 
	$idx=$1;
#    foreach my $idx (keys %names) {
	$line =~ s/$idx(\:|\s)/$names{$idx}$1/;
    }
    #$line =~ s/x\:/\:/g;
    print OUT "$line";
}
close OUT;
