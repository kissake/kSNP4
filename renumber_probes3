#!/usr/bin/perl
#v3.0

# renumber_probes.pl $probe_file.[0-9][0-9][0-9]
no warnings 'deprecated';

my @probe_files=@ARGV;
my %probeID =();
my $count_probe_seqs=0;

open OUT,">probes";
print OUT "SEQ_ID\tPROBE_ID\tPROBE_SEQUENCE\n";
open TABLE,">probes.table";
print TABLE "LOCUS_ID\tSTRAND\tALLELE\tPROBE_ID\tPROBE_SEQUENCE\n";

my $target="seq_id";

my $count=0;

foreach my $file (@probe_files) {
    my $previous=-1;
    print "$file\n";
    my $table_file=$file.".table";
    open IN,"$table_file";
    while (my $line=<IN>) {
	chomp $line;
	
	if ($line !~ /PROBE_ID/) {
	    my ($locus_ID,$strand,$allele,$pID,$seq)=split/\s+/,$line;
	    
	    $seq =~ s/[^atcgATCG]//g; # remove non-ATCG characters
	    
	    if ( length($seq)>25 )  {
		
		if ($locus_ID =~ /(.*)_(\d+)$/) {
		    $target=$1;
		    $num=$2;
		    
		    if ($previous ne $num) {
			$count++;
			$previous = $num;
		    }
		    $locus_ID =$target."_".$count;
		}
		if (!defined $probeID{$seq}) {
		    $count_probe_seqs++;
		    $probeID{$seq}=$count_probe_seqs;
		    print OUT "$target\t$probeID{$seq}\t$seq\n";
		}
		print TABLE "$locus_ID\t$strand\t$allele\t$probeID{$seq}\t$seq\n";
	    } # if ( length($seq)>25 )
	} # if ($line !~ /PROBE_ID/) {
    } #  while (my $line=<IN>) {
    close IN;
} # foreach $file 
