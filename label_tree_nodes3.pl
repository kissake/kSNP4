#!/usr/bin/perl
#v3.1

# This inputs a "tree", parses it with library, prints all of the node IDs,
# and then changes the non-leaf nodes to be identified by a number indicating
# the sequential order they were parsed into, then displays / outputs the
# complete tree (after the lists of the IDs?) to stdout. - JN

# ~/script/label_tree_nodes.pl outtree > outtree_nodes_labeled
no warnings 'deprecated';
use Bio::TreeIO;
use Bio::TreeIO::newick;

my $intree=$ARGV[0];

my $treeio = new Bio::TreeIO(-file   => "$intree",
                            -format => "newick");
my $tree = $treeio->next_tree;

@nodes=$tree->get_nodes();
my $count_nodes=0;
#print @nodes;
foreach $node (@nodes) {
    my $id=$node->id;
    #print "$id\n";
    if ($node->is_Leaf) {
	#Don't change the name
    } else {
    #if (!defined $id) {
	# Relabel with numbers. Delete any previous label.
	$count_nodes++;
	$node->id($count_nodes);
    }   
    
}

my  $treeio_nodes_named = new Bio::TreeIO( -format => "newick");

print $treeio_nodes_named->write_tree($tree),"\n";




