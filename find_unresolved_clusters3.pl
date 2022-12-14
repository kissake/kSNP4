#!/usr/bin/perl 
#v3.1

# find_unresolved_clusters.pl  NewickTreefile 

no warnings 'deprecated';
use Bio::TreeIO::newick;
use Net::FTP::Tiny;
use Bio::TreeIO;

my $treefile=$ARGV[0];

my $input = new Bio::TreeIO(-file   => "$treefile",
                            -format => "newick");
my $tree = $input->next_tree;

my @taxa = $tree->get_leaf_nodes;

my %strain=();

foreach my $leaf_node (@taxa) {
   my $id=$leaf_node->id;
   $strain{$id}=1;
}

#my $total_branch_length= $tree->total_branch_length;    

my %unresolved_cluster=();

foreach my $leaf_node1 (@taxa) {
    foreach my $leaf_node2 (@taxa) {
	my $dist=$tree->distance(-nodes => [$leaf_node1,$leaf_node2]);
	if ($leaf_node1 ne $leaf_node2 && $dist==0) {
	    my $strain1=$leaf_node1->id;
	    my $strain2=$leaf_node2->id;
            $unresolved_cluster{$strain1}{$strain2}=1;
        }
    }
}


#open OUT, ">unresolved_clusters" or warn "Cannot open unresolved_clusters file for writing: $!\n";
print "#Unresolved clusters\n";
my $cluster_count=0;
my %done=();
foreach my $nodeid1 (sort keys %unresolved_cluster) {
    if (!defined $done{$nodeid1}) {
        $cluster_count++;
        print  "$cluster_count\t$nodeid1\n";
        $done{$nodeid1}=1;
        foreach my $nodeid2 (sort keys %{$unresolved_cluster{$nodeid1}}) {
            print  "$cluster_count\t$nodeid2\n";
            $done{$nodeid2}=1;
        }
    }
}
print   "#Uniquely resolved\n";
foreach my $id1 (sort keys %strain){
    if (! %{$unresolved_cluster{$id1}}) {
        $cluster_count++;
        print  "$cluster_count\t$id1\n";
    }
}




