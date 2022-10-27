#!/usr/bin/perl
#v3.1

# core_SNPs.pl SNPs_all_labelLoci fileName2genomeName $min_fraction_with_locus
# creates core_SNPs and nonCore_SNPs and SNPs_in_majority.#
# core_SNPs has the SNP loci that are present in all the input genomes
# nonCore_SNPs has the rest of the loci, that are present in a subset of genomes
# SNPs_in_majority.$min_fraction_with_locus  has the SNPs present in at least this fraction of genomes: $min_fraction_with_locus
no warnings 'deprecated';

my $snp_file=$ARGV[0];
my $idfile=$ARGV[1];
my $min_fraction_with_locus=$ARGV[2];

open IN, "$idfile";
my $num_ids=0;
while (my $line = <IN>) {
    chomp($line);
    if ($line =~ /\w+/) {
	$num_ids++;
    }
}
close IN;

#print "$num_ids\n";
open IN, "$snp_file" or die "Cannot open $snp_file: $!.\n";
my $count=0;
my $locus_count=0;

my %by_id=();
my %snps=(); 
my $count_core=0;
my $count_noncore=0;
my $count_majority=0;
my $out="SNPs_in_majority".$min_fraction_with_locus;
open CORE,">core_SNPs" or warn "Cannot open core_SNPs for writing: $!\n";
open NONCORE,">nonCore_SNPs" or warn "Cannot open nonCore_SNPs for writing: $!\n";
open MAJORITY,">$out" or warn "Cannot open $out for writing: $!\n";

while (my $line = <IN>) {
    chomp($line);
   # print "$line\n";
    if ($line !~ /\w+/) {
	$locus_count++;
	$count=0;
	if ($locus_count >1) {
	    # print SNP into either core or nonCore file
	    my $num_ids_with_locus=scalar keys %snps;
	    if ($num_ids_with_locus >= $num_ids) {
		$count_core++;
		print CORE "\n";
		print MAJORITY  "\n";
		foreach my $id (sort {$a cmp $b} keys %snps) {
		    foreach my $c (sort {$a <=> $b} keys %{$snps{$id}}) {
			print CORE "$snps{$id}{$c}\n";
			print MAJORITY "$snps{$id}{$c}\n";
		    }
		}
	    } elsif  ($num_ids_with_locus >= $num_ids*$min_fraction_with_locus && $num_ids_with_locus < $num_ids) {
		$count_noncore++;
		$count_majority++;
		print MAJORITY  "\n";
		print NONCORE "\n";
		foreach my $id (sort {$a cmp $b} keys %snps) {
		    foreach my $c (sort {$a <=> $b} keys %{$snps{$id}}) {
			print MAJORITY "$snps{$id}{$c}\n";
			print NONCORE "$snps{$id}{$c}\n";
		    }
		}
	    } else {
		$count_noncore++;
		print NONCORE "\n";
		foreach my $id (sort {$a cmp $b} keys %snps) {
		    foreach my $c (sort {$a <=> $b} keys %{$snps{$id}}) {
			print NONCORE "$snps{$id}{$c}\n";
		    }
		}
	    }
	}
	# reinitialize %snp after printing
	%snps=(); 
    } elsif ($line =~ /\d+\t.*\t.*\t(\d+(?:\sF|\sR)?|x)\t(\S+)/) { 
	my $id=$2;
	$count++;
	$snps{$id}{$count}=$line;  # keep by count in case locus occurs multiple positions in a genome
    }
    
}

# Add last one
my $num_ids_with_locus=scalar keys %snps;
if ($num_ids_with_locus >= $num_ids) {
    $count_core++;
    print CORE "\n";
    print MAJORITY  "\n";
    foreach my $id (sort {$a cmp $b} keys %snps) {
	foreach my $c (sort {$a <=> $b} keys %{$snps{$id}}) {
	    print CORE "$snps{$id}{$c}\n";
	    print MAJORITY "$snps{$id}{$c}\n";
	}
    }
} elsif  ($num_ids_with_locus >= $num_ids*$min_fraction_with_locus && $num_ids_with_locus < $num_ids) {
    $count_noncore++;
    $count_majority++;
    print MAJORITY  "\n";
    print NONCORE "\n";
    foreach my $id (sort {$a cmp $b} keys %snps) {
	foreach my $c (sort {$a <=> $b} keys %{$snps{$id}}) {
	    print MAJORITY "$snps{$id}{$c}\n";
	    print NONCORE "$snps{$id}{$c}\n";
	}
    }
} else {
    $count_noncore++;
    print NONCORE "\n";
    foreach my $id (sort {$a cmp $b} keys %snps) {
	foreach my $c (sort {$a <=> $b} keys %{$snps{$id}}) {
	    print NONCORE "$snps{$id}{$c}\n";
	}
    }
}

$count_majority += $count_core;

close IN;
close CORE;
close NONCORE;

open COUNT,">COUNT_coreSNPs" or warn "Cannot open COUNT_coreSNPs: $!.\n";
print COUNT "Number core SNPs: $count_core\n";
print  COUNT "Number non-core SNPs: $count_noncore\n";
print  COUNT "Number SNPs in at least a fraction $min_fraction_with_locus of genomes: $count_majority\n";

close COUNT or warn "Cannot close COUNT_coreSNPs: $!.\n";