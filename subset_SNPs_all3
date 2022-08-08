#!/usr/bin/perl
#v3.0

#  subset_SNPs_all SNPs_all 25
#  subset_SNPs_all SNPs_all
# k is optional

no warnings 'deprecated';

my $infile=$ARGV[0]; # probably SNPs_all
my $k=5;
my @bases=qw(A C G T);
if ($ARGV[1]) {
    $k=$ARGV[1];
}

if ($k > 5) {
    my $count=0;
    foreach $w (@bases) {
	foreach $x (@bases) {
	    foreach $y (@bases) {
		foreach $m (@bases) {
		    foreach $n (@bases) {
			#foreach $t qw(A C G T) {
			#foreach $s qw(A C G T) {
			my $z=$w.$x.$y.$m.$n;
			$zarray[$count]=$z;
			$count++;
			#}
			#}
		    }
		}
	    }
	}
    }
} elsif ($k == 5) {
    my $count=0;
    foreach $w (@bases) {
	foreach $x (@bases) {
	    foreach $y (@bases) {
		foreach $m (@bases) {
		    #foreach $t qw(A C G T) {
		    #foreach $s qw(A C G T) {
		    my $z=$w.$x.$y.$m;
		    $zarray[$count]=$z;
		    $count++;
		    #}
		    #}
		}
	    }
	}
    }
} elsif ($k == 4) {
    my $count=0;
    foreach $w (@bases) {
	foreach $x (@bases) {
	    foreach $y (@bases) {
		my $z=$w.$x.$y;
		$zarray[$count]=$z;
		$count++;
	    }
	}
    }
} elsif ($k < 4) {
    my $count=0;
    foreach $w (@bases) {
	my $z=$w;
	$zarray[$count]=$z;
	$count++;
    }
}

$count=0;
open IN,"$infile";
my $out=$zarray[$count].".mers.SNPs_all";
open OUT,">$out";
while (my $line=<IN>){

    if ($line =~ /\S+/) {
	if ($line =~ /^\d+\t$zarray[$count]/) {
	    print OUT "$line";
	} else {
	    close OUT;
	    until ($line =~ /^\d+\t$zarray[$count]/) {
		$count++;
		$out=$zarray[$count].".mers.SNPs_all";
		open OUT,">$out";
		close OUT;
		
	    }
	    
	    $out=$zarray[$count].".mers.SNPs_all";
	    open OUT,">$out";
	    print OUT "$line";
	    
	}
    } # if ($line =~ /\S+/) {
}
close OUT;
