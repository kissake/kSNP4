#!/usr/bin/perl
#v3.0

#  labelTree_AlleleCount.pl tree_nodeLabel Node_SNP_counts  tree_AlleleCounts tree_nodeAllelecounts
#  labelTree_AlleleCount.pl tree_nodeLabel Node_SNP_counts  tree_AlleleCounts tree_nodeAllelecounts 1
# defaults to not number nodes (last arg 0), just show SNP counts. 
no warnings 'deprecated';
use Bio::TreeIO;

my $intree=$ARGV[0];
my $countFile=$ARGV[1];
my $new_treefile=$ARGV[2];
my $new_treefile_sametips=$ARGV[3];
my $number_nodes=0;
if ($ARGV[4] > 0 ) {
    $number_nodes=1;
}

my %nh=();
open IN,"$countFile";
my @lines=<IN>;
close IN;
chomp @lines;
foreach my $line (@lines) {
    if ($line =~ /node\:\s(.*)\sNumberTargets\:\s(\d+)\sNumber\S+\:\s(\d+)/) {
	$node=$1;
	$count_alleles=$3;
	$nh{$node}=$count_alleles;
    }
}

#foreach my $nodeID (keys %nh) {
#    print "$nodeID\n";
#}
#sleep(3);

my $treeio = new Bio::TreeIO(-file   => "$intree",
                            -format => "newick");
my $tree = $treeio->next_tree;

@nodes=$tree->get_nodes();

#print @nodes;
foreach $node (@nodes) {
    my $oldid=$node->id;
    if (!defined $nh{$oldid}) {
	$nh{$oldid}=0;
    }
    my $newid="";

    if ($node->is_Leaf) {
	$newid=$oldid."_".$nh{$oldid};
    } elsif ( $number_nodes== 0) {
	$newid=$nh{$oldid};
    } else {
	$newid=$oldid."_".$nh{$oldid};
    }

   # print "old id: $oldid\n";
   # print "new id: $newid\n";
    $node->id($newid);
    
}

my  $treeio_nodes_named = new Bio::TreeIO( -file => ">$new_treefile",
                                       -format => "newick");


$treeio_nodes_named->write_tree($tree);

$treeio = new Bio::TreeIO(-file   => "$intree",
                            -format => "newick");
$tree = $treeio->next_tree;
@nodes2=$tree->get_nodes();

foreach $node (@nodes2) {
    my $oldid=$node->id;
    if (!defined $nh{$oldid}) {
	$nh{$oldid}=0;
    }
    my $newid="";

    if ($node->is_Leaf) {
	$newid=$oldid;
    } elsif ( $number_nodes== 0) {
	$newid=$nh{$oldid};
    } else {
	$newid=$oldid."_".$nh{$oldid};
    }

   # print "old id: $oldid\n";
   # print "new id: $newid\n";
    $node->id($newid);
    
}



$treeio_nodes_named = new Bio::TreeIO( -file => ">$new_treefile_sametips",
                                       -format => "newick");


$treeio_nodes_named->write_tree($tree);

