#!/usr/bin/perl
#v3.1

# SNPs2nodes-new.pl SNPs_all_labelLoci nodes.perlhash FastTree Node_SNP_counts 
no warnings 'deprecated';

use Bio::TreeIO;
use Bio::TreeIO::newick;
use Storable;


my %clades=();

my $all_snps_file = $ARGV[0];
# A pre-computed hash that contains all possible clade structures
# for the given nodes.
my $NodeHashFile = $ARGV[1];
my $treefile=$ARGV[2];
my $NodeSigCounts=$ARGV[3];

# Able to get clade info from legacy text format file or reading in
# a pre-generated hash table.  Which is faster?  Would it be faster
# or otherwise better to put generator in the same script (no need
# to serialize the data at all) - JN
# 
# Determine file format from filename.  Potential for error? - JN
if ($NodeHashFile =~ /perlhash/) {
    %clades=%{ retrieve($NodeHashFile) };
} else {
    # it's a list, not stored perlhash, so parse the list
    open IN,"$NodeHashFile";
    my $node;
    my $root;

    while (my $line=<IN>) {
	chomp $line;
	if ($line =~ /root\:\s(\S+)\tnode\:\s(.*)/) {
	    $root=$1;
	    $node = $2;

	    print "node: $node\n";
	} elsif ($line =~ /\S/) {
	    push(@{$clades{$root}{$node}},$line);
	    print "id in node: $line\n";
	}
    }
}

print "Reading SNPs from file saved from prior run.\n";
my %fors = (); #hash storing snps
my %locus=();
my %all_ids=();
&read_snps_from_SNPs_all_file;
my $num_ids= scalar(keys %all_ids);

my %idsInNode=();
my %group=();
my %gl=();
my %core=();


# Traverse all nodes in %clades hash table (n squared) to create
# a %idsInNode hash table to permit referencing a node by the combo
# of the root and the list of nodes that leafed under it
# specifically.  Perhaps this should be loaded too, as it was
# created with the %clades hash. - JN

# @{$clades{$nodeID}}=($id1,$id2,....)
foreach $root (keys %clades) {
    foreach $node (keys %{$clades{$root}}) {
	my @node_ids_sorted=sort ( @{$clades{$root}{$node}} );
	my $idstring=join(",",@node_ids_sorted);
	$idsInNode{$idstring}{$root}=$node;
    }
}

my $total_SNPs=scalar (keys %fors);
print  "Number SNPs : $total_SNPs\n";

$date = `date`; print "$date\n";
    
my %num_map_to_node=();
my $snp_chars;

# Process the sequences in the order they occur in the SNPs_all file
# (does order matter?  If not, and if this SNPs_all file can get big,
# then the sorting may have a performance impact. - JN
foreach my $sequence ( sort {$locus{$a} <=> $locus{$b} } keys %fors) {  
    my @all_variants=();
    my %all_vars=();
    # Process the sequence in order by the ID of the genome source.
    foreach my $id (sort keys %{$fors{$sequence}}) {
	
	my %already_there=();
	# Cycle through all of the possible positions where this
	# sequence occurs in this genome (not sorted?), and for
	# each valid SNP, add the current genome to the list of
	# genomes where that SNP occurs for this sequence.
	foreach my $position ( keys %{$fors{$sequence}{$id}}) {
	    $snp_chars=$fors{$sequence}{$id}{$position};
	    if ($snp_chars eq "") {
		$snp_chars = "-";
	    }
	    #print "snp_chars: $snp_chars\n";
	    if (!defined $already_there{$snp_chars} ) { 
		# This ensures this genome is only added to the list
		# once for each SNP.
		push(@{$all_vars{$snp_chars}},$id);
		$already_there{$snp_chars} = 1;
	    }
	}

    }

    # make hash of Core SNPs
    if ( $num_ids > scalar(keys %{$fors{$sequence}}) ) {
	# This sequence doesn't occur in all genome sources.
	$core{$sequence}=0;
    } else {
	# This sequence DOES occur in all genome sources.
	$core{$sequence}=1;
    }

    # For each possible genome source
    foreach my $root (keys %clades) {
	my $yes_node=0;
	foreach my $snp_chars (keys %all_vars) {
	    # Recommend making join...sort...genomes into a subroutine - JN
	    my $idstring=join(",",sort(@{$all_vars{$snp_chars}}));
	    # Create an entry in the IDstring => sequence hash table %group,
	    $group{$idstring}{$sequence}=1;
	    # and create one in the sequence => IDstring hash table %gl
	    $gl{$sequence}{$idstring}=1;
	    if (defined $idsInNode{$idstring}{$root}) { # Not homoplastic
		# If this sequence + SNP occurs in all the leaves under any
		# node when this currently selected node was the root, this
		# is a good match, and a point in favor of this node being
		# the root. - JN
		$yes_node=1;
	    } 
	}
	if ($yes_node == 1) {
	    $num_map_to_node{$root}++;
	} # else {
	   # if ( $save_homoplasy  == 1) {
	#	foreach my $snp_chars (keys %all_vars) {
	#	    my $idstring=join(",",sort(@{$all_vars{$snp_chars}}));
	#	}
	 #   }
	#} # if ($yes_node == 0) {
    } # foreach my $root (keys %clades) {
} # foreach my $sequence 


my $max_num_node_loci=0;
my $best_root=-99;
print "Finding best root node\n";
foreach my $root (keys %clades) {

    if ($num_map_to_node{$root} > $max_num_node_loci) {
	$max_num_node_loci = $num_map_to_node{$root}  ; 
	$best_root=$root;
    }
}
print  "Number SNPs : $total_SNPs\n";

my $num_homoplastic_loci= $total_SNPs - $max_num_node_loci;

print "best root: $best_root\n";
print "number of homoplastic loci: $num_homoplastic_loci\n";
my $treeio = new Bio::TreeIO(-file   => "$treefile",
                            -format => "newick");
my $tree = $treeio->next_tree;

my $best_root_object=$tree->find_node(-id => "$best_root");
$tree->reroot($best_root_object);

my $new_treefile=$treefile.".rerooted";

my $treeio_rerooted = new Bio::TreeIO( -file => ">$new_treefile",
				       -format => "newick");

$treeio_rerooted->write_tree($tree);

print "Done writing tree\n";
`date`;

#my %in_node=();  # commented out July 24, 2012, see comment below

open COUNT,">COUNT_Homoplastic_SNPs";
print COUNT "Number_Homoplastic_SNPs: $num_homoplastic_loci\n";
close COUNT;

my %Clusters=();
print "Printing to $NodeSigCounts file\n";
open NODES,">$NodeSigCounts";

    
foreach my $nodeID ( sort keys %{$clades{$best_root}}) {
    
    my @ids=sort @{$clades{$best_root}{$nodeID}};
    my $count_ids=scalar @ids;
    my $idstring=join(",",@ids);
    my  $count_SNPs=0;
    if ($count_ids == 1) {
	$Clusters{$idstring}="Leaf.Node.".$nodeID;
    } else {
	$Clusters{$idstring}="Internal.Node.".$nodeID;
    }
    $count_SNPs=scalar keys(%{$group{$idstring}});
   # if ($count_SNPs > 0) {   # commented out July 24, 2012, see comment below
   #	foreach my $sequence (sort {$locus{$a} <=> $locus{$b} } keys %{$group{$idstring}} ) {
   #	    $in_node{$sequence}=1;
   #	}
   # }
    print NODES "node: $nodeID\tNumberTargets: $count_ids\tNumberSNPs: $count_SNPs\n";
    foreach my $id (@{$clades{$best_root}{$nodeID}}) {
	print NODES "$id\n";
    }
    print NODES "\n";
}

close NODES;
print "Done printing to $NodeSigCounts file\n";
`date`;


my $number_nodes= scalar (keys %{$clades{$best_root}});


#if ($save_homoplasy  == 1) {
    my %hp_ids=();
    my %hp_count=();

    # Start the unique numbering of homoplasy groups at the number of nodes (why?
    # Possibly so that we can number the nodes individually as their own groups /
    # have a number that is unique between nodes and groups?)
    my $homoplasy_group_count=$number_nodes;
    
    foreach my $idstring (sort {scalar(keys %{$group{$b}}) <=> scalar(keys %{$group{$a}}) } keys %group) {
	if (!defined $idsInNode{$idstring}{$best_root}) { # Homoplastic
	    my  $count_SNPs=0;
	    foreach my $sequence ( sort keys %{$group{$idstring}}) {
	#	if (!defined $in_node{$sequence}) { #no other alleles at this locus map to a node
	# Commented if statement out because I decided to make a homoplastic group for the alleles that do not map to a node
	# even if some of the alleles do map to a node. Previously, if any allele mapped to a node, no homoplastic group
	# was created for any of the alleles, even those that did not map. This is a change made on July 24, 2012
		    $count_SNPs++;
	#	}
	    }
	    
	  
	    if ($count_SNPs > 0) {
		# Uniquely number the homoplasy groups.
		$homoplasy_group_count++;
		$hp_ids{$homoplasy_group_count}=$idstring;
		$hp_count{$homoplasy_group_count}=$count_SNPs;

		$Clusters{$idstring}="Group.".$homoplasy_group_count;
	    

	    } #  if ($count_SNPs > 0) {
	} # 	if (!defined $idsInNode{$idstring}{$best_root}) {
    } # foreach my $idstring
    my $Homoplasy_groups="Homoplasy_groups";
    open GROUPS,">$Homoplasy_groups";

    # Sort the groups based on the value reported as NumberSNPs ($count_SNPs below)
    foreach my $homoplasy_group_count (sort { $hp_count{$b} <=> $hp_count{$a} } keys %hp_count) {   
	my $idstring=$hp_ids{$homoplasy_group_count};
	my $group="Group.".$homoplasy_group_count;
	my @these_ids=split/,/,$idstring;
	my $count_ids=scalar(@these_ids);
	$count_SNPs=$hp_count{$homoplasy_group_count};
	print GROUPS "Group: $group\tNumberTargets: $count_ids\tNumberSNPs: $count_SNPs\n";
	foreach my $id (@these_ids) {
	    print GROUPS "$id\n";
	}
	print GROUPS "\n";
    }
    close GROUPS;
    print "Done printing to $Homoplasy_groups\n";
    
#} # if ($save_homoplasy  == 1) {



my $outfile="ClusterInfo";
open OUT1,">$outfile";
print OUT1 "LocusNumber\tContextSeq\tCore\tClusters\n";

foreach my $sequence ( sort {$locus{$a} <=> $locus{$b} } keys %fors) {  
    my @clusterPrint=();
    foreach my $idstring (sort keys %{$gl{$sequence}}) {
	push @clusterPrint,$Clusters{$idstring};
    }
   # my $clustPrint=join(",",@clusterPrint);
    printf OUT1 "$locus{$sequence}\t$sequence\t$core{$sequence}\t@clusterPrint\n";
#    foreach my $id (sort keys %{$fors{$sequence}}) {
#	foreach my $position (sort { $a cmp $b} keys %{$fors{$sequence}{$id}}) {
#	    printf OUT1 "$locus{$sequence}\t$sequence\t$fors{$sequence}{$id}{$position}\t$position\t$id\t$core{$sequence}\t@clusterPrint\n";
#	}
#    }
#    print OUT1 "\n";
}


close OUT1;


$date = `date`; print "\n$date";
print "Finished SNPs2nodes-new\n";



##########################################################
########################## Subroutines ##############

# This subroutine is only used once. - JN
# Some of the structures generated by this subroutine (fors, locus)
# are generated in other contexts also.  Perhaps this can be a
# shared subroutine / module. - JN
#
# how do we cope when position data isn't available (set to "x")?
# It looks like that means that each k-mer can only occur once for
# that ID, or the snp_char will be overwritten and one or more SNP
# will be lost...?  - JN
# 
sub read_snps_from_SNPs_all_file {
    if (-e $all_snps_file) {
	
	open ALL_SNPS, "$all_snps_file" or die "Cannot open all $all_snps_file: $!\n";

	while (my $line = <ALL_SNPS>){
	    
	    #if ($line =~ /(.*)\t(.*)\t(\d+\s(?:F|R)|x)\t(.*)/) {  
	    if ($line =~ /(\d+)\t(.*)\t(.*)\t(\d+(?:\sF|\sR)?|x)\t(\S+)/ ) {
		#print "$1\t$2\t$3\t$4\t$5\n";
		my $locus_num=$1;
		my $primer = $2;
		my $snp_char = $3;
		my $position = $4;
		my $id= $5;
		chomp($id);
		
		$fors{$primer}{$id}{$position} = $snp_char;
		$locus{$primer}=$locus_num;
		$all_ids{$id}=1;

		#print "$primer\t$fors{$primer}{$id}{$position}\t$position\t$id\n";
	    }
	}
	
    }
} # sub read_snps_from_SNPs_all_file {




# This subroutine not used. - JN
# This subroutine also exists in at least one script (or similar).  Perhaps it should
# also be in a module - JN
###################
sub revcomp {
    my $seq = shift;
    my $comp = $seq;
    
    $comp =~ tr/ATGCatgc/TACGtacg/;
    my $revcomp_seq = reverse($comp);
    return $revcomp_seq;
}
#####################

