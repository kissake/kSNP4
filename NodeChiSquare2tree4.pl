#!/usr/bin/perl
#v3.1

=begin
significant SNPs are those with chi-square probabilities < p, where p is given in the second argument
Usage: NodeChiSquare -p 0.0001  -t ML
=cut
use warnings;
use strict;
use Statistics::Distributions;
use Cwd;
no warnings 'deprecated';
use Bio::TreeIO;
use Bio::TreeIO::newick;
use LWP::Simple;
use Getopt::Long; 

my $max_probability=0;
my $intree="";
my $help="";

GetOptions ( "f:s" =>\$intree,
             "p:f" =>\$max_probability,
	     "h" =>\$help);

my $nodeSNPcountfile = $intree;
$nodeSNPcountfile =~ s/tree_AlleleCounts(.*).tre/Node_SNP_counts$1/;

if ($intree eq "tree_AlleleCounts.parsimony.tre"){
    $nodeSNPcountfile = "Node_SNP_counts.SNPs_all.parsimony"
}

if ( not -e $intree or not -e $nodeSNPcountfile ){
    $help = 1;
    print "ERROR: Missing one or both of '$intree' or '$nodeSNPcountfile'.  Both are required.\n\n";
}

if ( $max_probability == 0){
    $help = 1;
    print "ERROR: Missing maximum probability option -p.\n\n"
}

if ($help) {
    print "Usage: NodeChiSquare2tree -p 0.001 -f tree_AlleleCounts.SNPs_all.parsimony.tre\n";
    print "\n";
    print " -f\tTree file to use as input (required) (assumes there is a corresponding Node_SNP_counts file) \n";
    print " -p\tMaximum Chi Sq probability for assigning a SNP to a node. (required)\n";
    print "\n";
    exit;
}

my $outtree=$intree;
$outtree=~s/tree_AlleleCounts/tree_ChiSqAlleleCounts/;

my $thisDate = localtime;
my $beginTime = time;
my $endTime = 0;
my $elapsedTime = 0;
my @args = ();
my @temp = ();
my $line = '';
my %Seq = (); #key is strain ID value is SNP sequence string
my @sequences = (); #index is SNP ID each element holds a sequence
my $numSeq = 0;
my $lenSeq = 0; #length of each sequence in @sequences
my @expecteds = (); #col0 = SNP ID, col 1 = prop A, col 2 = prop C, col3 = prop G, col4 = prop G, col 5 = prop gap
my @nodes = (); #col 0 = node ID, col 0 = array of names in the bipartition of that node
my $numNodes = 0;
my $aref; #pointer to the array of names stored in col 1 of @nodes

#read the SNPs_all_matrix.fasta file
readSNP_matrix(\%Seq, \@sequences);
$lenSeq = length $sequences[0];

#put the proportions of each allele for each SNP into @expecteds
getExpecteds (\@sequences, \@expecteds);
#print "$expecteds[1][0]\t$expecteds[1][1]\t$expecteds[1][2]\t$expecteds[1][3]\t$expecteds[1][4]\t$expecteds[1][5]\n";
$numSeq = scalar @sequences;

#parse the Node_SMNP_counts.ML file and store in @ nodes. col0 = node ID, col 1 areference to an anonymous array
readNode_SNP_counts ($numSeq, \@nodes);
$numNodes = scalar @nodes;
print "There are $numNodes nodes.\n";

#sort @nodes by node ID
@nodes = sort {$a-> [0] <=> $b-> [0]} @nodes;

#calculate and print the chi-square probabilities for significant SNPs at each node
open (OUTFILE, '>NodeChiSquares.txt') or die "Can't open NodeChiSquares.txt for writing. $!";
calculateChiSquares (\%Seq, \@nodes, $numSeq, $numNodes, $lenSeq, \@expecteds);

close OUTFILE;

&relabel_tree('NodeChiSquares.txt');

$endTime = time;
$elapsedTime = $endTime - $beginTime;
print "This run took $elapsedTime seconds.\n";
#############################################################################
                #Subroutines
#############################################################################
#called by Main
sub readSNP_matrix
{
my $Seq = $_[0]; #pointer to %Seq in which key is strain ID value is SNP sequence string
my $sequences = $_[1]; #pointer to @sequences
my $strainID = '';
my $line = '';
my @temp = ();



open(INFILE, 'SNPs_all_matrix.fasta') or die "Can't open SNPs_all_matrix.fasta for reading. $!";

while ($line = <INFILE>) {
	chomp $line;
	if ($line =~ /^>/) {
		@temp = split />/, $line;
		$strainID = $temp[1];
		$line = <INFILE>;
		chomp $line;
		push (@$sequences, $line);
		$$Seq{$strainID} = $line;
		}
	}

close INFILE;
}
#############################################################################
#Called by Main
sub getExpecteds
{
my $sequences = $_[0]; #pointer to @sequences
my $expected = $_[1]; #pointer to @expecteds
my $SNP = 0; #position in the SNP sequences
my $numSNPs = 0;
my $numSeq = 0;
my $numA = 0;
my $numC = 0;
my $numG = 0;
my $numT = 0;
my $numgap = 0;
my $propA = 0;
my $propC = 0;
my $propG = 0;
my $propT = 0;
my $propgap = 0;
my $theChar = '';

$numSeq = scalar @$sequences;
$numSNPs = length $$sequences[0];
#print "There are $numSeq sequences and $numSNPs SNPs in each sequence.\n";

for (my $i=0; $i< $numSNPs; $i++) {#
	for(my $j = 0; $j < $numSeq; $j++) {
		$theChar = substr($$sequences[$j], $i,1);
		if ($theChar eq 'A') {
			$numA++;
			}
			if ($theChar eq 'C') {
				$numC++;
				}		
				if ($theChar eq 'G') {
					$numG++;
					}	
					if ($theChar eq 'T') {
						$numT++;
						}		
						if ($theChar eq '-') {
							$numgap++;
							}
		}
	$propA = $numA/$numSeq;
	$propC = $numC/$numSeq;
	$propG = $numG/$numSeq;
	$propT = $numT/$numSeq;
	$propgap = $numgap/$numSeq;
	#print "$i\t$propA\t$propC\t$propG\t$propT\t$propgap\n";
	push @expecteds,[$i, $propA, $propC, $propG, $propT, $propgap];
	$numA = 0;
	$numC = 0;
	$numG = 0;
	$numT = 0;
	$numgap = 0;
	}

}
#############################################################################
#called by Main
sub readNode_SNP_counts
{
my $numSeq = $_[0];
my $nodes = $_[1]; #pointer to @nodes
my @temp = ();
my @temp2 = ();
my $line = '';
my $numTargets = 0;
my $nodeID = 0;
my @bipartition = (); #holds names of taxa in this bipartition
#my $aref; #anonymous array pointer




open (INFILE, $nodeSNPcountfile) or die "Can't open $nodeSNPcountfile for reading. $!";

while ($line = <INFILE>) {
	chomp $line;
	if ($line =~/^node:/) {
		@temp = split/\t/, $line;
		@temp2 = split / /, $temp[0];
		$nodeID = $temp2[1];
		@temp2 = split / /, $temp[1];
		$numTargets = $temp2[1];
		if ($numTargets >1 and $numTargets < $numSeq-1) {
			for(my $i = 0; $i < $numTargets; $i++) {
				$line = <INFILE>;
				chomp $line;
				push @bipartition, $line;
				}			
			#$aref = \@bipartition; #aref is a pointer to @ bipartition
			push @$nodes, [$nodeID,[@bipartition]];
			@bipartition = ();
			}
		}
	}
}

#############################################################################
#called by Main
sub calculateChiSquares
{
my $seq = $_[0]; #pointer to %seq
my $nodes = $_[1];
my $numSeq = $_[2];
my $numNodes = $_[3];
my $lenSeq = $_[4];
my $expecteds = $_[5]; #pointer to @expecteds col0 = SNP ID, col 1 = prop A, col 2 = prop C, col3 = prop G, col4 = prop G, col 5 = prop gap
my $theNode = 0;
my $aref; #pointer to the array that is stored in @nodes col 1
my $numStrains = 0; #number of strains in a node
my $theChar = '';
my $obsA = 0;
my $obsC = 0;
my $obsG = 0;
my $obsT = 0;
my $obsgap = 0;
my $expectedA = 0;
my $expectedC = 0;
my $expectedG = 0;
my $expectedT = 0;
my $expectedgap = 0;
my $chiA = 0;
my $chiC = 0;
my $chiG = 0;
my $chiT = 0;
my $chigap = 0;
my $chiSquare = 0;
my $probability = 0;
my $df = -1; #start at -1 because df = number of classes - 1
my @significantSNPs = (); #col0 = SNP ID, col 1 = probability
my $numSignificantSNPs = 0;


for (my $i = 0; $i < $numNodes; $i++) {#for each node i
	$theNode = $$nodes[$i][0];
	$aref = $$nodes[$i][1];
	$numStrains = scalar @$aref;
	#print "$theNode\t$numStrains\n";
	for(my $j = 0; $j < $lenSeq; $j++) {#for each site j in the sequence
		for( my $k = 0; $k < $numStrains; $k++) { #for each strain
			$theChar = substr($$seq{$$aref[$k]},$j,1);
			if ($theChar eq 'A') {
				$obsA++;
				}
				if ($theChar eq 'C') {
					$obsC++;
					}		
					if ($theChar eq 'G') {
						$obsG++;
						}	
						if ($theChar eq 'T') {
							$obsT++;
							}		
							if ($theChar eq '-') {
								$obsgap++;
								}					
			}
		#get the expecteds for this site
		$expectedA = $$expecteds[$j][1] * $numStrains;
		$expectedC = $$expecteds[$j][2] * $numStrains;
		$expectedG = $$expecteds[$j][3] * $numStrains;
		$expectedT = $$expecteds[$j][4] * $numStrains;
		$expectedgap = $$expecteds[$j][5] * $numStrains;
		
		#get the chi values for this site
		if( $expectedA > 0) {
			$chiA = (($obsA - $expectedA)**2)/$expectedA;
			$df++;
			}
		if( $expectedC > 0) {
			$chiC = (($obsC - $expectedC)**2)/$expectedC;
			$df++;
			}
		if( $expectedG > 0) {
			$chiG = (($obsG - $expectedG)**2)/$expectedG;
			$df++;
			}
		if( $expectedT > 0) {
			$chiT = (($obsT - $expectedT)**2)/$expectedT;
			$df++;
			}
		if( $expectedgap > 0) {
			$chigap = (($obsgap - $expectedgap)**2)/$expectedgap;
			$df++;
			}
		
		#calculate chiSquare
		$chiSquare = $chiA + $chiC + $chiG + $chiT + $chigap;
		
		#calculate the probability
		$probability = Statistics::Distributions::chisqrprob ($df,$chiSquare);
		
		#reset degrees of freedom to -1
		$df = 1;
		
		if ($probability<= $max_probability) {
			push @significantSNPs,[$j,$probability];
			}
		$obsA = 0;
		$obsC = 0;
		$obsG = 0;
		$obsT = 0;
		$obsgap = 0;
		$expectedA = 0;
		$expectedC = 0;
		$expectedG = 0;
		$expectedT = 0;
		$expectedgap = 0;
		$chiA = 0;
		$chiC = 0;
		$chiG = 0;
		$chiT = 0;

		}
	$numSignificantSNPs = scalar @significantSNPs;
	print OUTFILE "Node\t$theNode\tthere are\t$numSignificantSNPs\tsignificant SNPs at this node.\n";
	print OUTFILE "\tSNP\tProbability\n";
	for  (my $m = 0; $m < scalar @significantSNPs; $m++) {
		print OUTFILE "\t$significantSNPs[$m][0]\t$significantSNPs[$m][1]\n";
		}
	#reset the @significantSNPs aray to empty
	@significantSNPs = ();
	}
}

sub relabel_tree {
    my $countFile=shift;
    my %nh=();
    open IN,"$countFile";
    my @lines=<IN>;
    close IN;
    chomp @lines;
    foreach my $line (@lines) {
	if ($line =~ /Node\s+(\S+)\s+there\sare\s+(\d+)\s+significant/) { 
	    #print "$line\n";
	    my $node=$1;
	    my $count_alleles=$2;
	    #print "$node\n";
	    #print "$count_alleles\n";
	    $nh{$node}=$count_alleles;
	}
    }
    
    my $treeio = new Bio::TreeIO(-file   => "$intree",
			      -format => "newick");
    my $tree = $treeio->next_tree;
    my @nodes2=$tree->get_nodes();
    
    foreach my $node (@nodes2) {
	my $oldid=$node->id;
	if (!defined $nh{$oldid}) {
	    $nh{$oldid}=0;
	}
	my $newid="";
	
	if ($node->is_Leaf) {
	    $newid=$oldid;
	} else  {
	    $newid=$nh{$oldid};
	} 
	
	# print "old id: $oldid\n";
	# print "new id: $newid\n";
	$node->id($newid);
	
    }
    
    
    my $treeio_nodes_named = new Bio::TreeIO( -file => ">$outtree",
					   -format => "newick");
    $treeio_nodes_named->write_tree($tree);
    
}

#############################################################################

#############################################################################

#############################################################################

#############################################################################

#############################################################################

#############################################################################

#############################################################################
