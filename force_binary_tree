#!/usr/bin/perl
#v3.1
use Bio::TreeIO;
use Bio::TreeIO::newick;
use Bio::Phylo::Forest::Tree;
no warnings 'deprecated';

my $treefile=$ARGV[0];
my $outtree=$ARGV[1];

my $input = new Bio::TreeIO(-file   => "$treefile",
                            -format => "newick");
my $tree = $input->next_tree;
my $phylotree = Bio::Phylo::Forest::Tree->new_from_bioperl($tree);

my @nodes=$phylotree->get_nodes;
print "number nodes " , scalar @nodes,"\n";
$phylotree->resolve();
 @nodes=$phylotree->get_nodes;
print "number new nodes " , scalar @nodes,"\n";

$treeio = new Bio::TreeIO( -file => ">$outtree",
                                       -format => "newick");

$treeio->write_tree($phylotree);

exit;

