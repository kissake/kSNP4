#!/usr/bin/perl
#v3.0

# parse mummer output to find SNP position and strand in genome

# parse_mummer4kSNPv2.pl   $mummer_output > SNP.positions

no warnings 'deprecated';

my $mummer_out=$ARGV[0];

open IN,"$mummer_out";
my $id;
my $strand;

while (my $line=<IN>) {
    chomp $line;
    if ($line =~ /^>\s(\S+)/) {
	$id=$1;
	chomp($id);
	$strand="F";

	if ($line =~ /^>.*(Reverse)$/) {
	    $strand="R";
	}
	
    }
    
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

	if ($strand eq "F") {
	    $snp_pos=($hit_len-1)/2+$hit_start;
	} else {
	    $snp_pos=$hit_start - ($hit_len-1)/2;
	}
	if ($snp_locus_allele =~ /(\S+)_([ACTG])/) {
	    $snp_locus=$1;
	    $snp_allele=$2;
	}
	print "$snp_locus\t$snp_allele\t$snp_pos $strand\t$id\n";
    }
    
}
close IN;

