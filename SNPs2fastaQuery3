#!/usr/bin/perl
#v3.0

# Shea 10-19-04

# $Bin/SNPs2fastaQuery.pl $original_SNPs_all >! SNPloci.fasta
no warnings 'deprecated';

my $snp_file=$ARGV[0];

my %snps=(); # $snps{surrounding.sequence}{ snp character}=$line;

open IN, "$snp_file" or die "Cannot open $snp_file: $!.\n";

while (my $line = <IN>) {
    chomp($line);

    if ($line =~ /\d+\t(.*)\t(.*)\t(\d+(?:\sF|\sR)?|x)\t(.*)/) { 

	$count++;
	my $allele=$2;
	my $locus=$1;
	my ($s1,$s2)=split/\./,$locus;
	my $seq=$s1.$allele.$s2;
	$snps{$locus}{$allele}=$seq;

    }
}
close IN;

foreach my $locus  (sort keys %snps) {
    foreach my $allele ( sort keys %{$snps{$locus}}) {
	print ">$locus","_",$allele,"\n$snps{$locus}{$allele}\n";
    }
} 

