#!/usr/bin/perl
#v3.0


# parse_SNPs2VCF.pl SNPs_all_labelLoci $ref_genome $outfile
# choose the ref genome and name it as second argument exactly as it appears in the SNPs_all_labelLoci file

no warnings 'deprecated';

my %fors=();
my %locus=();
my $all_snps_file=$ARGV[0];
my $out=$ARGV[1];
my $ref_genome=$ARGV[2];

my $out_not_shown="VCF.SNPsNotinRef.".$ref_genome;

&read_snps_from_SNPs_all_file($all_snps_file);

$count_snps=scalar keys %fors;

print "\$count_snps: $count_snps\n";

my $chromosome=1;

my %count_genomes=();
foreach my $sequence ( sort keys %fors) {
    foreach my $id (sort keys %{$fors{$sequence}}) {
	$count_genomes{$id}=1;
	if ($ref_genome eq "") {
	    $ref_genome=$id;
	}
    }
}
my $genome_count=scalar keys %count_genomes;

open(OUT2, ">$out_not_shown") || die "Can't open [$out_not_shown] for writing $!\n";
open(OUT, ">$out") || die "Can't open [$out] for writing $!\n";
print OUT "##fileformat=VCFv4.0\n";
print OUT "##Reference genome=$ref_genome\n";
print OUT "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
print OUT "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n";

print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach my $id (sort keys %count_genomes) {
    print OUT "\t$id";
}
print OUT "\n";


my $count=0;
foreach my $sequence ( sort keys %fors) {
    my $onReverse=0;
    foreach my $position (sort {$a cmp $b} keys %{$fors{$sequence}{$ref_genome}}) {
	if ($position !~ /R/) {
	    $ref_allele{$sequence}=$fors{$sequence}{$ref_genome}{$position};
	} else {
	    $ref_allele{$sequence}=revcomp($fors{$sequence}{$ref_genome}{$position});
	    $onReverse=1;
	}
    }

    foreach my $id (sort keys %{$fors{$sequence}}) {
        foreach my $position (sort {$a cmp $b} keys %{$fors{$sequence}{$id}}) {
	    my $a;
	    if ($onReverse == 0) {
		$a=$fors{$sequence}{$id}{$position};
	    } else {
		$a=revcomp($fors{$sequence}{$id}{$position});
	    }

	    if ($a ne $ref_allele{$sequence} && $fors{$sequence}{$id}{$position} =~ /[ATCG]/) {
		$alt_allele{$sequence}{$a}++;
	    } 
	}
    }

    my @alt=();
    my @alt_freq=();
    foreach my $alt (sort keys %{$alt_allele{$sequence}}) {
	push @alt,$alt;
	my $af=sprintf("%.3f",$alt_allele{$sequence}{$alt}/$genome_count);
	push @alt_freq,$af;
    }
    my $alt_string=join(",",@alt);
    my $alt_freq_string=join(",",@alt_freq);
    my $num_samples=scalar keys %{$fors{$sequence}};
    my $strand="";
    foreach my $position (sort {$a cmp $b} keys %{$fors{$sequence}{$ref_genome}}) {
	if ($position =~ /F/) {
	    $strand="_F";
	} elsif ($position =~ /R/) {
	    $strand="_R";
	    $onReverse=1;
	}
	my $pos_noStrand="";
	if ($position =~ /(\d+)/) {
	    $pos_noStrand=$1;
	}

	#if ($position !~ /(\w)/  ) {
	 #   foreach my $id (sort keys %{$fors{$sequence}}) {
	#	foreach my $position (sort {$a cmp $b} keys %{$fors{$sequence}{$id}}) {
	#	    print "$sequence\t$fors{$sequence}{$id}{$position}\t$position\t$id\n";
	#	}
	  #  }
	#}	
    
	print OUT "$chromosome\t$pos_noStrand\t$sequence","$strand\t$ref_allele{$sequence}\t$alt_string\t.\t.\tNS=$num_samples;AF=$alt_freq_string\tGT";
	my %allele=();
	foreach my $id (sort keys %count_genomes) {
	    $allele{$id}=".";
	    foreach my $pos2 ( keys %{$fors{$sequence}{$id}}) {
		my $x="";
		if ($onReverse == 0) {
		    $x=$fors{$sequence}{$id}{$pos2};
		} else {
		    $x=revcomp($fors{$sequence}{$id}{$pos2});
		}
		if ($x  eq $ref_allele{$sequence}) {
		    $allele{$id}=0;
		} else {
		    for (my $a=0; $a < @alt ; $a++) {
			if ($x eq $alt[$a]) {
			    $allele{$id}=$a+1;
			}
		    }
		}
	    }
	    print OUT "\t$allele{$id}";
	} # foreach my $id (sort keys %count_genomes) {
	print OUT "\n";
    } #  foreach my $position 
    if (scalar(keys %{$fors{$sequence}{$ref_genome}}) < 1) {
	foreach my $id (sort keys %{$fors{$sequence}}) {
	    foreach my $position (sort {$a cmp $b} keys %{$fors{$sequence}{$id}}) {
		print OUT2 "$sequence\t$fors{$sequence}{$id}{$position}\t$position\t$id\n";
	    }
	}
	print OUT2  "\n";
    }
    
} # foreach my $sequence (sort keys %fors) 



close OUT;
close OUT2;



sub read_snps_from_SNPs_all_file {
    my  $all_snps_file=shift;
    open ALL_SNPS, "$all_snps_file" or die "Cannot open all $all_snps_file: $!\n";
    my $primer="X";
    while (my $line = <ALL_SNPS>){

        if ($line =~ /^(\d+)\t(.*)\t(.*)\t(\d+(?:\sF|\sR)?|x)\t(\S+)/ ) {
            my $locus_num=$1;
            $primer = $2;
            my $snp_char = $3;
            my $position = $4;
            my $name= $5;
            chomp($name);
            $fors{$primer}{$name}{$position} = $snp_char;   # global, but gets reset ### - JN - Beg to differ.  Does not get reset.
            $locus{$primer}=$locus_num;  # global
            #print "$primer\t$fors{$primer}{$id}{$position}\t$position\t$id\n";
        }
    }
    close ALL_SNPS or die "Cannot close $out_prefix.$suffix: $!.\n";
} # sub read_snps_from_SNPs_all_file {

sub revcomp {
    my ($s) = @_;
    $s =~ tr/wsatugcyrkmbdhvnATUGCYRKMBDHVN/WSTAACGRYMKVHDBNTAACGRYMKVHDBN/;
    $s = reverse $s;
    return $s;
}
