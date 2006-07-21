#!/usr/bin/perl -w
#
# This script replaces a copyright header in file with the one
# that is input.
#
use strict;
#
my $g_use_msg =
  "Use: update-copyright-header.pl new_copyright_file file_in file_out\n";
if( scalar(@ARGV) != 3 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
#
my $copyright_file    = shift;
my $file_in_name  = shift;
my $file_out_name = shift;
#
open COPYRIGHT_FILE, "<$copyright_file" || die "The file $copyright_file could not be opened!";
my @copyright_lines = <COPYRIGHT_FILE>;
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
my $found_first_header_statement = 0;
my $wrote_header = 0;
my @temp_file_out_array;
my @file_out_array;
#
foreach(@file_in_array) {
  my $line = $_;
  if(!$wrote_header) {
    if($line=~/\@HEADER/) {
      if(!$found_first_header_statement) {
        # This is the first HEADER tag so let's write the lines
        # we have read up to this point.
        push @file_out_array, @temp_file_out_array;
        $found_first_header_statement = 1;
      }
      else {
        # This is the second HEADER tag so lets insert the
        # new copyright header and just copy the rest of the file.
        push @file_out_array, @copyright_lines;
        $wrote_header = 1;
        next;
      }
    }
    if(!$found_first_header_statement) {
      # We have not found the first @HEADER tag so let's
      # copy the lines into the temp array until we do.
      push @temp_file_out_array, $line;
    }
  }
  else {
    # We have already written the header so we just need
    # to copy the rest of the lines in the file.
    push @file_out_array, $line;
  }
}

if(!$wrote_header) {
  # The header was never written so we will just
  # insert the copyright header at the top and
  # copy the rest
  push @file_out_array, @copyright_lines;
  push @file_out_array, @file_in_array;
  $wrote_header = 1;
}

if( $wrote_header || $file_out_name ne $file_in_name ) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
