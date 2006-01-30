#!/usr/bin/perl -w
#
# This perl script fixes up doxygen documenation formating
# from what rabartl was using to what is used more standardly
# in Trilinos.
# Here it is allowed for file_in and file_out to be the
# same file.
#
use strict;
#
my $g_use_msg =
  "Use: string-replace.pl find_string file_in file_out\n";
if( scalar(@ARGV) < 2 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
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
my $i;
foreach( $i = 0; $i < scalar(@file_in_array); ++$i ) {
  my $line = $file_in_array[$i];
  #print "line = $line";
  if( $line=~m[^\s*///$] ) {
    $did_replacement = 1;
    my $next_line = $file_in_array[$i+1];
    #print "next_line = $next_line";
    if( $next_line=~m[/\*\*] ) {
      $next_line=~s[/\*\*][/** \\brief];
      push @file_out_array, $next_line;
      ++$i;
    }
    else {
      $line =~ s[///][/** \\brief . */];
      push @file_out_array, $line;
    }
  }
  else {
    $did_replacement = 1 if $line =~ s[\\brief \\brief][\\brief];
    push @file_out_array, $line;
  }
  #print "line added = $file_out_array[-1]";
}
if($did_replacement || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
