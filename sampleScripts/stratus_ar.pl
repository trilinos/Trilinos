#!/usr/bin/perl
# This script was originally developed to use on DEC platforms that have
# a argument size limitation.  For templated codes, the arguments were so
# long that the archiver complained of "arg list too long".  Instead this
# script will parse the archiver line and deal with the cxx_repository, if
# necessary.  Otherwise, nothing will change with the archiving line.
#
# Author : Heidi Thornquist
# Date : May 19, 2004
#

$output = "";
$i = 0;
$idx = 0;
$numinputs = @ARGV;
$templates = 0;
#
# Check to see if there is any template instantiations for this archiving.
# [ $idx will hold the index to the first cxx_repository statement when this finishes ]
#
while($idx < $numinputs) {
  if (index($ARGV[$idx],"cxx_repository")==0) {
    $templates = 1;
    last;
  }
  $idx++;
}
if ($templates) {
  #
  # First get the archiver, options, and library name.
  #
  $ar_opts_lib = "";
  $j = 0;
  while( index($ARGV[$j],".o") < 0 ) {
    $ar_opts_lib = $ar_opts_lib . $ARGV[$j] . " ";
    $j++;
  }
  #
  # First archive the non-template object files made by the compiler.
  #
  while($i < $idx) {
    $output = $output . $ARGV[$i] . " ";
    $i++;
  }
  system "$output";    # Execute the archiving of the non-template object files.
  #
  # Reset the output line to have the archiving executable, options, and the library.
  #
  $file_pattern = "";
  @pattern_list = ("[a-m]*", "[A-M]*", "[n-z]*", "[N-Z]*", "_*");
  foreach $pattern (@pattern_list) {
    $file_pattern = "cxx_repository/" . $pattern; 
    @file_list = glob($file_pattern);
    if (@file_list > 1) {
      $output = $ar_opts_lib . " ";
      foreach $filename (@file_list) {
	$output = $output . $filename . " ";
      }
      system "$output";  # Execute the archiving of this set of template object files.
    }
  }
} else {
  #
  # Just pass all the arguments through, there are no template instantiations to deal with.
  #
  while ($i < $numinputs) {
    $output = $output . $ARGV[$i] . " ";
    $i++;
  }
  exec ("$output") || die "Cannot run reformatted archiving line";
}
# We should be done here!

