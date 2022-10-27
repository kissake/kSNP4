#!/usr/bin/perl
#v3.03

# SNPs_all_2_fasta_matrix.pl SNPs_all_labelLoci SNPs_all_matrix.fasta SNPs_all_matrix
no warnings 'deprecated';

my $snp_file=$ARGV[0];
my $out_fasta=$ARGV[1];
my $out_matrix=$ARGV[2];

open IN, "$snp_file" or die "Cannot open $snp_file: $!.\n";
my $count_of_positions=0;
my %by_id=();
my $previous_locus=-1;

while (my $line = <IN>) {
    chomp($line);
    # print "$line\n";
   
    if ($line =~ /(\d+)\t\S*\t(\S*)\t(\d+(?:\sF|\sR)?|x)\t(\S+)/) {  
	my $locus=$1;
	my $snp=$2;
	my $id=$4;
	if ($previous_locus != $locus) {
            $count_of_positions++;
            $previous_locus = $locus;
        }

	$by_id{$id}{$count_of_positions}=$snp;

    }
}
close IN;
my $seqlen =$count_of_positions; 
my $num_ids=scalar keys %by_id;

#print "$num_ids   $seqlen\n";
open FASTA,">$out_fasta" or die "Cannot open $out_fasta: $!.\n";
open MATRIX, ">$out_matrix" or die "Cannot open $out_matrix: $!.\n";
print MATRIX "$num_ids   $seqlen\n";

foreach my $id (sort {$a cmp $b} keys %by_id) {
    
    print FASTA ">$id\n";
    printf MATRIX  "%-90.90s",$id;
    for (my $i=1;$i<=$seqlen;$i++) {
	if (defined $by_id{$id}{$i} &&  $by_id{$id}{$i} =~ /[ATCG]/) {
	    $snp_fasta=$by_id{$id}{$i};
	    $snp_matrix=$by_id{$id}{$i};
	} else {
	   $snp_fasta = "-";
	   $snp_matrix = "N";
	}
    
	print FASTA "$snp_fasta";
	print MATRIX "$snp_matrix";
    }
    print FASTA "\n";
    print MATRIX  "\n";
}

close FASTA or die "Cannot close  $out_fasta: $!.\n";
close MATRIX or die "Cannot close  $out_matrix: $!.\n";
