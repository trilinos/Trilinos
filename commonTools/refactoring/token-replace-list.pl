#!/usr/bin/perl -w
#
# This perl script replaces a list of tokens
# on a token basis.  Here it is allowed for file_in and
# file_out to be the same file.
#
use strict;
#
my $g_use_msg =
  "Use: token-replace-list.pl token_file file_in file_out\n";
if( scalar(@ARGV) < 3 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
#
my $token_file    = shift;
my $file_in_name  = shift;
my $file_out_name = shift;
my $print_file_name = shift;
#
open TOKEN_NAMES_FILE, "<$token_file" || die "The file $token_file could not be opened!";
my @token_names = <TOKEN_NAMES_FILE>;
#
#print "file_in_name = $file_in_name\n";
if($file_in_name=~/CVS/) {
#  print "Do not replace in CVS\n";
  exit;
}
open FILE_IN, "<$file_in_name";
my @file_in_array = <FILE_IN>;
close FILE_IN;
#
#print $match_str . "\n";
#
my @file_out_array;
my $did_replacement = 0;
my $printed_file_name = 0;
foreach(@file_in_array) {
  my $line = $_;
  foreach(@token_names) {
	my ($find_token, $replace_token) = split(" ",$_);
	my $match_str = '([^a-zA-Z_]|^)' . $find_token . '([^a-zA-Z0-9_]|$)';
	$did_replacement = 1 if $line=~s/$match_str/$1$replace_token$2/g;
        if ( $did_replacement && $print_file_name && !$printed_file_name ) {
          print "Did token replacement in file : $file_in_name\n";
          $printed_file_name = 1;
        }
  }
  push @file_out_array, $line;
}
if($did_replacement || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
