#!/usr/bin/perl -w
#
# This perl script replaces the leading tabs in a
# text file with spaces.  Here it is allowd for
# file_in and file_out to be the same file.
#
use strict;
#
my $g_use_msg =
  "Use: replace-tabs-with-spaces spaces_per_tab file_in file_out\n";
if( scalar(@ARGV) != 3 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
#
my $spaces_per_tab = shift;
my $file_in_name   = shift;
my $file_out_name  = shift;
#
open FILE_IN, "<$file_in_name" || die "The file $file_in_name could not be opended for input\n";
my @file_in_array = <FILE_IN>;
close FILE_IN;
#
my $tab_char = '\t';
my @file_out_array;
my $found_leading_tab = 0;
foreach(@file_in_array) {
  my $len = length($_);
  my $i;
  my $line = "";
  my $end_of_leading_white_space = 0;
  for( $i = 0; $i < $len; ++$i ) {
    my $c = substr($_,$i,1);
    if( $c ne "\t" && $c ne " " ) {
      $end_of_leading_white_space = 1;
    }
    if( $end_of_leading_white_space || $c ne "\t" ) {
      $line .= $c;
      next;
    }
    $found_leading_tab = 1;
    for( my $j = 0; $j < $spaces_per_tab; ++$j ) {
      $line .= ' ';
    }
  }
  push @file_out_array, $line;
}
if($found_leading_tab || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
