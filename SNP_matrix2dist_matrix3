#!/usr/bin/perl 
#v3.0

# Shea Gardner Jan 5,2005
# Example: ~/script/SNP/SNP_matrix2dist_matrix.pl SNPs_all_matrix > NJ.dist.matrix

# echo "NJ.dist.matrix\ny" | neighbor -l 90
# 

no warnings 'deprecated';

open IN,"$ARGV[0]";
my %strain=();
while (my $line=<IN>) {
    chomp $line;
    if ($line =~ /(.{90})([\w-]+)/) {
	$id=$1;
	$seq=$2;
	$id =~ s/\s+//g;
	$strain{$id}=$seq;
    }
}

my %dist=();
foreach my $id1 (sort keys %strain){
    foreach my $id2 (sort keys %strain){
	$dist{$id1}{$id2}=0;
    }
}

foreach my $id1 (sort keys %strain){
    foreach my $id2 (sort keys %strain) {
	@seq1=split//,$strain{$id1};
	@seq2=split//,$strain{$id2};
	#print "len $id1: ", scalar @seq1,"\n";
	#print "len $id2: ", scalar @seq2,"\n";
	
	for (my $i=0; $i < @seq1; $i++) {
	    if (!defined $seq2[$i] ) {
		print "$seq2[$i]\n";
		print "$i\t$id1\t$id2\t$seq1[$i]\n";
	    }
	    if ($seq1[$i] ne $seq2[$i] && ($seq1[$i] eq "N" || $seq2[$i]  eq "N")  ) {
		$dist{$id1}{$id2} +=2;
	    } elsif ($seq1[$i] ne $seq2[$i] && $seq1[$i] ne "N" && $seq2[$i]  ne "N"  ) {
		$dist{$id1}{$id2} +=1;
	    } elsif ( $seq1[$i] eq "N" && $seq2[$i]  eq "N"  ) {
		$dist{$id1}{$id2} +=0.0;
	    }
	}
    }
}



my $num_seqs=keys %strain;
print " $num_seqs\n";
foreach my $id1 (sort keys %strain){
    #print "$id1 DIST ";
    my $id10char=$id1;
    $id10char =~ s/^(\S+)\s.*/$1/;
    $id10char=sprintf("%-90.90s",$id10char); 
    print "$id10char";
    	
    foreach my $id2 (sort keys %strain) {
	print "$dist{$id1}{$id2} ";
    }
    print "\n";
}




