#!/usr/bin/perl -w
#
# This perl script replaces a string with another string
# on a token basis.  Here it is allowed for file_in and
# file_out to be the same file.
#
use strict;
#
my $g_use_msg =
  "Use: token-replace.pl find_token replacement_token file_in file_out\n";
if( scalar(@ARGV) < 4 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
#
my $find_token         = shift;
my $replacement_token  = shift;
my $file_in_name       = shift;
my $file_out_name      = shift;
#
#print "file_in_name = $file_in_name\n";
if($file_in_name=~/CVS/) {
#  print "Do not replace in CVS\n";
  exit;
}
open FILE_IN, "<$file_in_name" || die "The file $file_in_name could not be opended for input\n";
my @file_in_array = <FILE_IN>;
close FILE_IN;
#
my $match_str = '([^\w\d_]|^)' . $find_token . '([^\w\d_]|$)';
#print $match_str . "\n";
#
my @file_out_array;
my $did_replacement = 0;
foreach(@file_in_array) {
  $did_replacement = 1 if $_=~s/$match_str/$1$replacement_token$2/g;
  push @file_out_array, $_;
}
if($did_replacement || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
