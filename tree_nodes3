#!/usr/bin/perl
#v3.1

#
# first label nodes: 
#$Bin/label_tree_nodes.pl tree.$t   > ! tree_nodeLabel.$t
# example:  $Bin/tree_nodes.pl tree_nodeLabel.$t   nodes.$t  

use Bio::TreeIO;
use Bio::TreeIO::newick;
use Storable;

no warnings 'deprecated';

my $treefile=$ARGV[0];
my $NodeFile=$ARGV[1];
my $outfile_unresolved_clusters=$ARGV[2];
my $NodeHashFile=$NodeFile.".perlhash";

my $input = new Bio::TreeIO(-file   => "$treefile",
                            -format => "newick");
my $tree = $input->next_tree;


my @nodes=$tree->get_nodes;

my $original_root = $tree->get_root_node;
#print "\$original_root: $original_root\n";
my @taxa = $tree->get_leaf_nodes;
my $count_leaves=scalar @taxa;

#my $half_node=$count_leaves*0.5;
#$half_node=int($half_node);
#$half_node = $tree->find_node(-id => "$half_node");

##print "\$half_node: $half_node \n";
#$quarter_node=$count_leaves*0.25;
#$quarter_node=int($quarter_node);
#$quarter_node = $tree->find_node(-id => "$quarter_node");
#$threequarter_node=$count_leaves*0.75;
#$threequarter_node=int($threequarter_node);
#$threequarter_node = $tree->find_node(-id => "$threequarter_node");

#@nodes2reroot=($original_root ,$half_node ,$quarter_node ,$threequarter_node );
@nodes2reroot=@nodes;  # @taxa;
my %clades_by_root=();
my %node_id_hash=();
my %clades=();

foreach my $n (@nodes2reroot) {
    #print "\n\n\n\n\nnode2reroot: ",$n->id,"\n";
    $tree->reroot($n);
    $root=$n->id;

    foreach my $node (@nodes) {
	#print   "\nnode: " , $node->id,"\n";
	#print $node->id, " get_all_Descendents\n";
	#print "Node number $count_nodes\n";

	if ($node->is_Leaf) {
	    $this_id=$node->id;
	    @{$clades{$node->id}}=($this_id);
	    @{$clades_by_root{$root}{$node->id}}=($this_id);
	} else {
	    for my $child ( $node->get_all_Descendents ) {
		if ($child->is_Leaf) {
		    push(@{$clades_by_root{$root}{$node->id}},$child->id);
		} 
	    }
	    @{$clades_by_root{$root}{$node->id}} = sort (@{$clades_by_root{$root}{$node->id}});
	    $idstring=join(",",sort (@{$clades_by_root{$root}{$node->id}}) );
	    $node_id_hash{$root}{$idstring}=$node->id;
	}
    }
}



open OUT,">$NodeFile";
foreach my $root (sort keys %clades_by_root) {
    foreach $node (sort { $a cmp $b  || $a <=> $b } keys %{$clades_by_root{$root}}) {  
	print OUT  "\nroot: $root\tnode: $node\n";
	foreach my $childID (@{$clades_by_root{$root}{$node}}) {
	    print OUT "$childID\n";
	}
    }
}
close OUT;
%clades=%clades_by_root;
store(\%clades,$NodeHashFile);

#print "number of roots: ",scalar(keys %clades_by_root),"\n";

print "number of nodes: ",scalar(keys %{$clades_by_root{$original_root->id}} ),"\n";


exit;





 
