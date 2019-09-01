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
 
    # $seq is never converted into canonical capitalized format, which could
    # impact this comparison?! TODO - JN
    my @array=sort($seq,$rc);
   # print "@array\n";
    $kmer{$array[0]} += $freq;

}
close IN;

# This 'sort' in the foreach is probably painful.  Faster in perl, or in shell?
# Also, if in shell is faster, there is an opportunity to reduce code.
foreach my $seq (sort keys %kmer) {
    print  "$seq $kmer{$seq}\n";
}


 
# Note, this 'reverse' step is non-trivial, results in all characters being
# capitalized, U's being lossily converted into A's, W, S, and N's not getting 
# transformed, and other values being swapped (e.g. A for T). - JN
sub revcomp {
    my ($s) = @_;
    $s =~ tr/wsatugcyrkmbdhvnATUGCYRKMBDHVN/WSTAACGRYMKVHDBNTAACGRYMKVHDBN/;
    $s = reverse $s;
    return $s;
}
