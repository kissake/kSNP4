#!/usr/bin/perl 
#v3.0

#delete_allele_conflicts.pl $kmer_counts 
# delete kmers from this genome in wihch the central base differs but the surrounding context is the same.

no warnings 'deprecated';

my $kmer_counts_file=$ARGV[0];

my $out=$kmer_counts_file.".conflictsDeleted";


my $seqLine=`head -1 $kmer_counts_file `;
my $seq;
my $k;
if ($seqLine=~ /(\w+)\s+\d+/) {
    $seq=$1;
    $k=length($seq);
}


my %kmers=();

open IN,"$ARGV[0]";
while (my $line=<IN>) {
    chomp $line;
   # print "line: $line\n";
    my ($seq,$freq)=split/\s+/,$line;
   # print "seq=$seq\n";
   # print "freq: $freq\n";
    my $s1=substr($seq,0,(($k-1)/2));
 
    my $s2=substr($seq,(($k-1)/2)+1);
    my $center=substr($seq,(($k-1)/2),1);
    my $j=$s1.".".$s2; # snp
    #print "$j\t$center\t$freq\n";
  
    $kmers{$j}{$center}=$freq;

}
close IN;

open OUT,">$out";
my %keep=();
foreach my $seq (sort keys %kmers ) {
    
    @center= keys(%{$kmers{$seq}});
    if ( scalar @center  == 1) {
	my ($s1,$s2)=split/\./,$seq;
	my $kseq= $s1.$center[0].$s2;
	$keep{$kseq}=1;
	#print OUT $s1,$center[0],$s2,"\t$kmers{$seq}{$center[0]}\n";
    }
}
foreach my $kseq (sort keys %keep) {
    print OUT "$kseq\n";
}

close OUT;
