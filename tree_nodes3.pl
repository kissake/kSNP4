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


#
# Reroot (create a new tree of the same nodes with the same connections, but 
# rooted in a specific node) using every node as a potential root node.
#
foreach my $n (@nodes2reroot) {
    #print "\n\n\n\n\nnode2reroot: ",$n->id,"\n";
    $tree->reroot($n);
    $root=$n->id;

    # In the newly re-rooted tree, for each node in the tree,
    foreach my $node (@nodes) {
	#print   "\nnode: " , $node->id,"\n";
	#print $node->id, " get_all_Descendents\n";
	#print "Node number $count_nodes\n";
	
	# Either the node is a leaf, or it has leaves...
	if ($node->is_Leaf) {
	    # If it is a leaf, then it is the sole member of its branch, so 
	    # the clades that branch off of it only include itself.
	    $this_id=$node->id;
	    @{$clades{$node->id}}=($this_id);
	    @{$clades_by_root{$root}{$node->id}}=($this_id);
	} else {
	    # If it has leaves, traverse the entire tree beneath it, adding only
	    # the ultimate leaves to the %clades_by_root hash.
	    for my $child ( $node->get_all_Descendents ) {
		if ($child->is_Leaf) {
		    push(@{$clades_by_root{$root}{$node->id}},$child->id);
		} 
	    }
	    # Then, sort the list of leaves (note that in the above case, a
	    # list of one element is already sorted) to create the ID string
	    @{$clades_by_root{$root}{$node->id}} = sort (@{$clades_by_root{$root}{$node->id}});
	    $idstring=join(",",sort (@{$clades_by_root{$root}{$node->id}}) );

	    # Then add that ID string to a new hash table to be able to create
	    # a backreference from the ID string for this root to the node that
	    # is associated with it.
	    $node_id_hash{$root}{$idstring}=$node->id;
	}
    }
}


#
# Output information about the structure of the nodes in the tree in the form
# of listing all possible combinations of two nodes, and listing all of the 
# leaf nodes of the second (in sort order) when the first is the root.
# 
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

#
# Then store the clades_by_root data (same as in nodes.parsimony) in a file in
# pre-computed form (no need to re-build the hash table, but you could from the
# data in the nodes.parsimony file).  That filename is nodes.parsimony.perlHash. 
%clades=%clades_by_root;
store(\%clades,$NodeHashFile);

#print "number of roots: ",scalar(keys %clades_by_root),"\n";

print "number of nodes: ",scalar(keys %{$clades_by_root{$original_root->id}} ),"\n";


exit;





 
