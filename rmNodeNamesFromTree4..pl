#!/usr/bin/perl
#v3.1

# rm_node_names_from_tree.pl tree_nodeLabel.$t.tre tree.$t.tre

no warnings 'deprecated';
use Bio::TreeIO;
use IO::Socket::SSL;
use Bio::TreeIO::newick;


my $intree=$ARGV[0];
my $new_treefile=$ARGV[1];

my $treeio = new Bio::TreeIO(-file   => "$intree",
                            -format => "newick");
my $tree = $treeio->next_tree;

@nodes=$tree->get_nodes();

#print @nodes;
foreach $node (@nodes) {
    my $oldid=$node->id;

    my $newid="";

    if ($node->is_Leaf) {
        $newid=$oldid;
    } else {
        $newid="";
    }

   # print "old id: $oldid\n";
   # print "new id: $newid\n";
    $node->id($newid);
    
}

my  $treeio_nodes_unnamed = new Bio::TreeIO( -file => ">$new_treefile",
                                       -format => "newick");


$treeio_nodes_unnamed->write_tree($tree);

