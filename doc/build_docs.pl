#!/usr/bin/perl

###############################################################################
# Trilinos/doc/build_docs.pl
#
# - You must run this script from this directoy!
# - Run any build_docs in any doc directory
# - Create html file with links to each set of documentation
#
###############################################################################

use strict;
use Getopt::Long;

# Read in commandline options
my $skipBuild = 0;

GetOptions(
  "skip-build!" => \$skipBuild
   );

# Figure out Trilinos directory
my $trilinosDir = `pwd`;
chomp($trilinosDir);
$trilinosDir =~ s/\/doc$//;

# Array for list of html/index.html paths
my @indexList = ();

# Run any build_docs script in any doc directory
# Fill array of html/index.html paths
use File::Find;
use File::Basename;
print "\n";
if(!$skipBuild) {
  print "\nGenerating package documentation ...\n";
  find (\&buildDocs, $trilinosDir);
}
else {
  print "\nSkipping generation of package documentation ...\n";
}
print "\nGenerating list of documentation directories ...\n";
find (\&fillIndexList, "$trilinosDir/packages");

#==============================================================================
# Generate links file

print "\nGenerating Trilinos PDF files ...\n";

# Open index.html for list of links to generated documentation
open(PAGE, ">index.html") or die "can't open index.html, died";

# Header information
print PAGE "<html><head><title>Trilinos Documentation</title></head><body>\n";

# Trilinos documentation section
print PAGE "<h3>Trilinos Documentation:</h3>\n";
print PAGE "<ul>\n";

# Variables for use in generating pdfs
my $output = "";
my $failed = 0;

# Trilinos Overview
print "Trilinos Overview...\n";
chdir("$trilinosDir/doc/OverviewSANDreport");
if(!$skipBuild) {
  $output = `make -f classicMakefile pdf 2>&1`;
  $failed = $?;
}
if (!$failed) { print PAGE "<li><a href=\"./OverviewSANDreport/TrilinosOverview.pdf\">Trilinos Overview (.pdf)</a></li>\n"; }

# Trilinos User Guide
print "Trilinos User Guide...\n";
chdir("$trilinosDir/doc/UserGuide");
if(!$skipBuild) {
  $output = `make -f classicMakefile pdf 2>&1`;
  $failed = $?;
}
if (!$failed) { print PAGE "<li><a href=\"./UserGuide/TrilinosUserGuide.pdf\">Trilinos User Guide (.pdf)</a></li>\n"; }

# Trilinos Tutorial
print "Trilinos Tutorial...\n";
chdir("$trilinosDir/packages/didasko/doc");
if (stat("tutorial.pdf")) { print PAGE "<li><a href=\"../packages/didasko/doc/tutorial.pdf\">Trilinos Tutorial (.pdf)</a></li>\n"; }

# Trilinos Epetra Performance Optimization Guide
print "Trilinos Overview...\n";
chdir("./packages/epetra/doc/PerfOptGuide");
if (stat("EpetraPerformanceGuide.pdf")) { print PAGE "<li><a href=\"./packages/epetra/doc/PerfOptGuide/EpetraPerformanceGuide.pdf\">Epetra Performance Guide (.pdf)</a></li>\n"; }

# Trilinos Developer Guide
print "Trilinos Developer Guide...\n";
chdir("$trilinosDir/doc/DevGuide");
if(!$skipBuild) {
  $output = `make -f classicMakefile pdf 2>&1`;
  $failed = $?;
}
if (!$failed) { print PAGE "<li><a href=\"./DevGuide/TrilinosDevGuide.pdf\">Trilinos Developers Guide (.pdf)</a></li>\n"; }

# Package documentation section
print PAGE "</ul>\n";
print PAGE "<h3>Trilinos Package Documentation:</h3>\n";
print PAGE "<ul>\n";

# Print sorted list of html/index.html links
foreach my $item (sort { lc($a) cmp lc($b) } @indexList) {
    print PAGE $item;
}

# Footer information
print PAGE "</ul>\n";
print PAGE "</body></html>\n";

# Close index.html
close(PAGE) or die "can't close index.html, died";

# Inform user how to access docs
print "\nFinished generating Trilinos documentation.\n";
print "Open ./index.html in your browser.\n\n";

###############################################################################
# Subroutines

# Run any build_docs script in any doc directory
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

# Fill array of html/index.html paths
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
        if($shortDir ne "") {
          push(@indexList, "<li><a href=\"../packages/$shortDir/doc/html/index.html\">$shortDir</a></li>\n");
        }
    
    }    
    
} # fillIndexList
