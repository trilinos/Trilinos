#!/usr/bin/perl -w
#
# This perl script replaces a list of strings.
# Here it is allowed for file_in and
# file_out to be the same file.
#

use strict;

my $g_use_msg =
  "Use: string-replace-list.pl string_file file_in file_out\n";
if( scalar(@ARGV) != 3 ) {
  print STDERR $g_use_msg;
  exit(-1);
}

my $string_file   = shift;
my $file_in_name  = shift;
my $file_out_name = shift;

open STRING_NAMES_FILE, "<$string_file" || die "The file $string_file could not be opened!";
my @string_names_lines = <STRING_NAMES_FILE>;
my @string_names;
foreach(@string_names_lines) {
  my $line = $_;
  my @string_names_line_array = split(' ', $line);
  foreach(@string_names_line_array) {
    chomp($_);
    push @string_names, $_
  }
}

#print "file_in_name = $file_in_name\n";
if($file_in_name=~/CVS/) {
#  print "Do not replace in CVS\n";
  exit;
}
open FILE_IN, "<$file_in_name";
my @file_in_array = <FILE_IN>;
close FILE_IN;

#print $match_str . "\n";

my @file_out_array;
my $did_replacement = 0;
foreach(@file_in_array) {
  my $line = $_;
  #print "\n";
  #print "----------------------------------------------------\n";
  #print "\n";
  #print $line;
  my $i;
  for( $i=0; $i <scalar(@string_names); $i += 2 ) {
    my ($find_string, $replace_string) = ($string_names[$i], $string_names[$i+1]);
    defined($find_string) || die $!;
    defined($replace_string) || die $!;
    chomp($find_string);
    #print "find string = ${find_string}\n";
    chomp($replace_string);
    #print "replace string = ${replace_string}\n";
    #print "string_names line = $_";
    $did_replacement = 1 if $line=~s/$find_string/$replace_string/g;
    #print $line;
  }
  push @file_out_array, $line;
}
if($did_replacement || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
