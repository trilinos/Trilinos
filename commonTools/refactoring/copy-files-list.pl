#!/usr/bin/perl -w
#
# This perl script reads a set of files and copies them from
# on directory to another and renames them given a list of names.
#
use strict;
#
my $g_use_msg =
  "Use: copy-files-list.pl old_dir new_dir old_to_new_file_list\n";
if( scalar(@ARGV) != 3 ) {
  print STDERR $g_use_msg;
  exit(1);
}
#
my $old_dir              = shift;
my $new_dir              = shift;
my $old_to_new_file_list = shift;
#print "old_dir = \'$old_dir\'\n";
#print "new_dir = \'$new_dir\'\n";
#print "old_to_new_file_list = \'$old_to_new_file_list\'\n";
#
open FILE, "<$old_to_new_file_list" || die "The file $old_to_new_file_list could not be opened!";
my @old_to_new_file_list = <FILE>;
close FILE;
#
foreach(@old_to_new_file_list) {
	my ($old_file_name, $new_file_name) = split;
	#print "old_file_name = \'$old_file_name\'\n";
	#print "new_file_name = \'$new_file_name\'\n";
	my $cmnd = "cp ${old_dir}/${old_file_name} ${new_dir}/${new_file_name}";
	#print "${cmnd}\n";
	system($cmnd)==0 || die "Error, could not executed the command \'$cmnd\': $?";
}
