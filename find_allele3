#!/usr/bin/perl
#v3.0


# find_allele.pl SNP_loci.fasta genome_kmer_file genome_name
no warnings 'deprecated';

my $snp_fasta=$ARGV[0];
my $genome_kmer_file=$ARGV[1];
my $genome_name=$ARGV[2];
my %kmers=();

open IN,"$genome_kmer_file";
while (my $line=<IN>) {
    chomp $line;
    $kmers{$line}=1;
    #print "$line\n";
}
close IN;

open IN, "$snp_fasta";
while (my $line=<IN>) {
    chomp $line;
    if ($line =~ /^>(.*)_([ACTG])/) {
	$locus=$1;
	$allele=$2;
	$seq=<IN>;
	chomp $seq;
	#print "locus: $locus\tallele: $allele\tseq: $seq\n";
	#$present=`fgrep -c  $seq $genome_kmer_file`;
	#print "present: $present\n";
	#if ($present ==1) {
	if (defined $kmers{$seq}) { #
	    print "$locus\t$allele\t$seq\t$genome_name\n";
	}
    }
}
