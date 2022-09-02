#!/usr/bin/perl
#v3.0

my $f=$ARGV[0];
my $this_dir=$ARGV[1];
#print "orig $f\n";
my $newfile="";
chomp $f;
chomp $this_dir;

if ($f =~ /^[\/|~]/ ) {
    # it already has path
    $newfile=$f;
} else {
    $newfile = $this_dir."/".$f;
}
print "$newfile\n";