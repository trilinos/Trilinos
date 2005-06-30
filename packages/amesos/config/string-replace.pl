#!/usr/bin/perl -w
#
# This perl script replaces a string with another string.
# Here it is allowd for file_in and file_out to be the
# same file.
#
use strict;
#
my $g_use_msg =
  "Use: string-replace.pl find_string replacement_string file_in file_out\n";
if( scalar(@ARGV) < 4 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
#
my $find_string        = shift;
my $replacement_string = shift;
my $file_in_name       = shift;
my $file_out_name      = shift;
#
#
if($file_in_name=~/CVS/) {
#  print "Do not replace in CVS\n";
  exit;
}
#
open FILE_IN, "<$file_in_name" || die "The file $file_in_name could not be opended for input\n";
my @file_in_array = <FILE_IN>;
close FILE_IN;
#
my @file_out_array;
my $did_replacement = 0;
foreach(@file_in_array) {
  #print $_;
  $did_replacement = 1 if $_=~s/$find_string/$replacement_string/g;
  #print $_;
  push @file_out_array, $_;
}
if($did_replacement || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
