#!/usr/bin/perl

################################################################################
# Trilinos/doc/build_docs.pl
#
# - Run any build_docs in doc directory
# - Create html file with links to each set of documentation
#
################################################################################

use strict;

my $trilinosDir = `pwd`;
chomp($trilinosDir);
$trilinosDir =~ s/\/doc$//;

my @indexList = ();

use File::Find;
use File::Basename;
find (\&buildDocs, $trilinosDir);
find (\&fillIndexList, $trilinosDir);

open(PAGE, ">index.html") or die "can't open index.html, died";

print PAGE "<html><head><title>Trilinos Documentation</title></head><body>\n";
print PAGE "<h3>Trilinos Documentation</h3>\n";
print PAGE "<h4>Package Documentation:</h4>\n";
print PAGE "<ul>\n";

foreach my $item (sort { lc($a) cmp lc($b) } @indexList) {
    print PAGE $item;
}

print PAGE "</ul>\n";
print PAGE "</body></html>\n";

close(PAGE) or die "can't close index.html, died";

print "\nFinished generating Trilinos documentation.\n";
print "Open ./index.html in your browser.\n\n";

sub buildDocs {
    
    use File::Basename;
    my $absFile = $File::Find::name;
    my $file = basename($absFile);
    my $dir = dirname($absFile);
    
    # Is this file in a 'doc' directory?
    # Is this file named 'build_docs'?
    if ($dir =~ m/doc$/ && $file =~ m/build_docs$/) {
    
        my $shortDir = $absFile;
        $shortDir =~ m/.*\/packages\/(.*)\/doc.*/;
        $shortDir = $1;
        print "$shortDir...\n";
        
        my $output = `$absFile 2>&1`;
        my $failed = $?;
    
    }    
    
} # buildDocs

sub fillIndexList {
    
    use File::Basename;
    my $absFile = $File::Find::name;
    my $file = basename($absFile);
    my $dir = dirname($absFile);
    
    # Is this file in an 'html' directory?
    # Is this file named 'index.html'?
    if ($dir =~ m/html$/ && $file =~ m/^index.html$/) {
        
        my $shortDir = $absFile;
        $shortDir =~ m/.*\/packages\/(.*)\/doc.*/;
        $shortDir = $1;
        
        push(@indexList, "<li><a href=\"$absFile\">$shortDir</a></li>\n");
    
    }    
    
} # fillIndexList