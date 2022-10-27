#!/usr/bin/perl
#v3.0

#Example: rc_kmer_freqs.pl kmer_counts_file_from_sa > kmer_counts_correcting4rc

no warnings 'deprecated';

my %kmer=();
open IN,"$ARGV[0]";

while (my $line=<IN>) {
    chomp($line);
    #$progress++;
   # if ($progress % 10000 == 0) { 
       # $date=`date`;
        #print  "$date  on line $progress\n";
    #}
    ($seq,$freq)=split/\s+/,$line;

    $rc=revcomp($seq);
    # only store the canonical kmer, one that appears first in sorted order, but count the frequency for occurence in both directions.
 
    my @array=sort($seq,$rc);
   # print "@array\n";
    $kmer{$array[0]} += $freq;

}
close IN;

foreach my $seq (sort keys %kmer) {
    print  "$seq $kmer{$seq}\n";
}


 
sub revcomp {
    my ($s) = @_;
    $s =~ tr/wsatugcyrkmbdhvnATUGCYRKMBDHVN/WSTAACGRYMKVHDBNTAACGRYMKVHDBN/;
    $s = reverse $s;
    return $s;
}

