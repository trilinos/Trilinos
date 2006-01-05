#!/usr/bin/perl -w
#
# This perl script is used to setup a build directory for
# code that is external to Trilinos but uses the hard work
# of the Trilinos autconf/automake system to define how
# code is compiled, packaged into libraries and linked
# into executables.
#
#
#
# Note, this script must be maintained to be current for
# the Teuchos makefile.
#
use strict;
use Cwd;
use File::Path;
#
# Pares the command-line
#
if( !defined(@ARGV) || scalar(@ARGV) < 3 ) {
  die
    "Usage: external-configure.pl TRILNOS_BUILD_DIR EXTERNAL_SRC_DIR PACKAGE1 PACAKGE2 ...\n".
    "  TRILNOS_BUILD_DIR : The path to the base Trilinos build directory\n".
    "  EXTERNAL_SRC_DIR  : The path to the base directory for PACKAGE1, PACAKGE2 ...\n"
    ;
}
my $trilinos_build_dir = make_abs_path(cwd(),shift);
my $external_src_dir = make_abs_path(cwd(),shift);
my @external_packages = @ARGV;
# ToDo: Validate that $trilinos_build_dir is an absolute path!
print "  trilinos_build_dir = $trilinos_build_dir\n";
print "  external_src_dir = $external_src_dir\n";
print "  external_packages = [", join(",",@external_packages), "]\n";
#
# Get the source directory for Trilinos
#
my $srcdir_line = `grep \'^srcdir.*=\' $trilinos_build_dir/Makefile`;
#print "srcdir_line: $srcdir_line";
my $trilinos_src_dir = make_abs_path($trilinos_build_dir,(split(" = ",$srcdir_line))[1]);
chomp $trilinos_src_dir;
# ToDo: If above is a relative path then we must append it to 
print "  trilinos_src_dir = $trilinos_src_dir\n";
#
# Create a base-level Makefile that will build the other projects
#
open BASE_MAKEFILE, ">Makefile";
print BASE_MAKEFILE ":\n";
foreach(@external_packages) {
  print BASE_MAKEFILE "\tcd ${_}; make\n";
}
print BASE_MAKEFILE "all :\n";
foreach(@external_packages) {
  print BASE_MAKEFILE "\tcd ${_}; make all\n";
}
print BASE_MAKEFILE "clean :\n";
foreach(@external_packages) {
  print BASE_MAKEFILE "\tcd ${_}; make clean\n";
}
print BASE_MAKEFILE "clean-obj :\n";
foreach(@external_packages) {
  print BASE_MAKEFILE "\tcd ${_}; make clean-obj\n";
}
print BASE_MAKEFILE "clean-lib :\n";
foreach(@external_packages) {
  print BASE_MAKEFILE "\tcd ${_}; make clean-lib\n";
}
print BASE_MAKEFILE "clean-exe :\n";
foreach(@external_packages) {
  print BASE_MAKEFILE "\tcd ${_}; make clean-exe\n";
}
#
# Create the external subpackage directories and their makefiles
#
foreach(@external_packages) {
  my $external_package = $_;
  print "  Setting up \"$external_package\" ...\n";
  (mkpath($external_package,1,0777) || die $!) if !(-e $external_package);
  run_cmnd("cp ${external_src_dir}/${external_package}/Makefile.in ${external_package}/Makefile",1);
  run_cmnd("${trilinos_src_dir}/commonTools/refactoring/token-replace.pl _TRILINOS_BUILD_DIR ${trilinos_build_dir} ${external_package}/Makefile ${external_package}/Makefile",1);
  run_cmnd("${trilinos_src_dir}/commonTools/refactoring/token-replace.pl _TRILINOS_SRC_DIR ${trilinos_src_dir} ${external_package}/Makefile ${external_package}/Makefile",1);
  run_cmnd("${trilinos_src_dir}/commonTools/refactoring/token-replace.pl _BASE_SRC_DIR ${external_src_dir}/${external_package} ${external_package}/Makefile ${external_package}/Makefile",1);
}

sub run_cmnd {
  my $cmnd = shift;
  my $stop_on_fail = shift;
  #print "$cmnd\n";
  if(system($cmnd) != 0 && $stop_on_fail) { die; }
}

sub make_abs_path {
  my $base_path = shift;
  my $path_in = shift;
  return "$base_path/$path_in" if( $path_in =~ /^\./ );
  return $path_in;
}
