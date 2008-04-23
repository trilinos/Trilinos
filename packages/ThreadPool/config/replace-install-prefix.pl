#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#
# This script is called to do a set of text replacements for installing
# a Mafile.export.package file so that external clients can use it.
#
# Read in commandline arguments
#
my $exec_prefix = "";           # [required] Abs path to base installation directory (i.e. --prefix=??? option passed to configure)
my $my_export_makefile = "";    # [required] Name only of installed Makefile.export.package file
my $my_top_srcdir = "";         # [required] Abs path to this package's top source directory
my $my_incl_dirs = "";          # [required] Abs path to this package's include directories
my $my_lib_dirs = "";           # [optional] Abs path to this package's library directories (if any exist)
my $dep_package_builddirs = ""; # [optional] Abs paths to other directly dependent framework package build directories (if any exist)
GetOptions(
  "exec-prefix=s"                   => \$exec_prefix,
  "my-export-makefile=s"            => \$my_export_makefile,
  "my-abs-top-srcdir=s"             => \$my_top_srcdir,
  "my-abs-incl-dirs=s"              => \$my_incl_dirs,
  "my-abs-lib-dirs=s"               => \$my_lib_dirs,
  "dep-package-abs-builddirs=s"     => \$dep_package_builddirs
  );
#
# Validate commandline arguments
#
scalar(@ARGV) == 0 || die;
$exec_prefix ne "" || die;
$my_export_makefile ne "" || die;
$my_top_srcdir ne "" || die;
$my_incl_dirs ne "" || die;
#
# Interpret commandline arguments
#
$exec_prefix = remove_rel_paths($exec_prefix);
my @my_incl_dirs = split(":",$my_incl_dirs);
my @my_lib_dirs = split(":",$my_lib_dirs);
my @dep_export_package_builddirs = split(":",$dep_package_builddirs);
#
# Do the replacements
#
my $my_abs_export_makefile = "${exec_prefix}/include/${my_export_makefile}";

my $cmnd_base = "${my_top_srcdir}/config/token-replace.pl ";
#
foreach(@dep_export_package_builddirs) {
  if($_ ne "") {
    run_cmnd($cmnd_base . "${_} ${exec_prefix}/include ${my_abs_export_makefile} ${my_abs_export_makefile}");
  }
}
#
foreach(@my_incl_dirs) {
  if($_ ne "") {
    run_cmnd($cmnd_base . "-I${_} -I${exec_prefix}/include ${my_abs_export_makefile} ${my_abs_export_makefile}");
  }
}
#
foreach(@my_lib_dirs) {
  if($_ ne "") {
    run_cmnd($cmnd_base . "-L${_} -L${exec_prefix}/lib ${my_abs_export_makefile} ${my_abs_export_makefile}");
  }
}
#
run_cmnd($cmnd_base . "${my_top_srcdir}/config ${exec_prefix}/include ${my_abs_export_makefile} ${my_abs_export_makefile}");
#
# Subroutines
#
sub remove_rel_paths {
	my $entry_in = shift;
	if ($entry_in=~/-L\.\./) {
		return $entry_in;
	}
	my @paths = split("/",$entry_in);
	my @new_paths;
	foreach( @paths ) {
		if( !($_=~/\.\./) ) {
			push @new_paths, $_;
		}
		else {
			pop @new_paths
		}
	}
	return join("/",@new_paths);
}
sub run_cmnd {
  my $cmnd = shift;
  #print "\n", $cmnd, "\n";
  system($cmnd)==0 || die;
}
