#!/usr/bin/perl
#v3.0

# parse mummer output to find SNP position and strand in genome

# parse_mummer4kSNPv2.pl   $mummer_output > SNP.positions

no warnings 'deprecated';

# Input is rows of 3 integers, with occasional '>'-prefixed lines for context? -JN
my $mummer_out=$ARGV[0];

# Output: tab separated list of records:
# <locus> <allele> <position><space><strand> <ID>
#
open IN,"$mummer_out";
my $id;
my $strand;

while (my $line=<IN>) {
    chomp $line;
    if ($line =~ /^>\s(\S+)/) {
	$id=$1;
	chomp($id);
	# Default to forward direction? - JN
	$strand="F";

	if ($line =~ /^>.*(Reverse)$/) {
            # If the word 'Reverse' is in the description, adjust to reverse 
	    # direction?  - JN
	    $strand="R";
	}
	
    }
    
    # Lines are:
    # <locus>_<allele><arbitrary non-whitespace> <unknown number> <hit_start> <hit_len>
    if ($line =~ /(\S+)\s+\d+\s+(\d+)\s+(\d+)/) {
	my $snp_locus_allele=$1;
	my $hit_start=$2;
	my $hit_len=$3;
	my $snp_pos;
	my $snp_locus;
	my $snp_allele;
	#print "\$hit_start: $hit_start  $strand\n";
	#print "\$hit_len: $hit_len\n";
	#print "$id\n";

	# If the strand is reversed, must travel backwards to reach the position
	# of the unique char for the k-mer - JN
	if ($strand eq "F") {
	    $snp_pos=($hit_len-1)/2+$hit_start;
	} else {
	    $snp_pos=$hit_start - ($hit_len-1)/2;
	}
	# No provision is made for this if statement to fail, in which case
	# snp_locus and snp_allele are null
	if ($snp_locus_allele =~ /(\S+)_([ACTG])/) {
	    $snp_locus=$1;
	    $snp_allele=$2;
	}
	print "$snp_locus\t$snp_allele\t$snp_pos $strand\t$id\n";
    }
    
}
close IN;

