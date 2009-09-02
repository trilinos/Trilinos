#!/usr/bin/perl -w
#
# This perl script removes the ^M from dos text files.
# Here it is allowd for file_in and file_out to be the
# same file.
#
use strict;
#
my $g_use_msg =
  "Use: fix-dos-text-file.pl file_in file_out\n";
if( scalar(@ARGV) != 2 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
#
my $file_in_name  = shift;
my $file_out_name = shift;
#
open FILE_IN, "<$file_in_name" || die "The file $file_in_name could not be opended for input\n";
my @file_in_array = <FILE_IN>;
close FILE_IN;
#
my $bad_char = 13;
my @file_out_array;
my $found_bad_char = 0;
foreach(@file_in_array) {
  my $len = length($_);
  my $i;
  my $line = "";
  for( $i = 0; $i < $len; ++$i ) {
	my $c = substr($_,$i,1);
	my $ascii_c = ord($c);
	if( $ascii_c == $bad_char ) {
	  $found_bad_char = 1;
	}
	else {
	  $line .= $c;
	}
  }
  push @file_out_array, $line;
}
if($found_bad_char || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
