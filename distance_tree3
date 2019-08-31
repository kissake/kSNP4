#!/usr/bin/perl
#v3.0

no warnings 'deprecated';

use Bio::Tree::DistanceFactory;
use Bio::TreeIO;
use Bio::Matrix::IO;

my $dfactory= Bio::Tree::DistanceFactory->new(-method=>'NJ');
my $treeout=Bio::TreeIO->new(-format=>'newick');


my $parser= Bio::Matrix::IO->new(-format => 'phylip',
				 -file => 'NJ.dist.matrix');

my $mat = $parser->next_matrix;

my $tree = $dfactory->make_tree($mat);

$treeout->write_tree($tree);
