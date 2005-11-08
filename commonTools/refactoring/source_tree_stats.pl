#!/usr/bin/perl -w
#
# Must be run from the base directory

use FileHandle;
use Cwd;

my $basedir = Cwd::cwd();

my ($num_hpp_files,  $num_hpp_lines, $num_hpp_semicolons) = get_files_lines_semicolons_counts("hpp");
my ($num_cpp_files, $num_cpp_lines, $num_cpp_semicolons)  = get_files_lines_semicolons_counts("cpp");
my ($num_H_files,   $num_H_lines, $num_H_semicolons)    = get_files_lines_semicolons_counts("H");
my ($num_C_files,   $num_C_lines, $num_C_semicolons)    = get_files_lines_semicolons_counts("C");
my ($num_h_files,   $num_h_lines, $num_h_semicolons)    = get_files_lines_semicolons_counts("h");
my ($num_c_files,   $num_c_lines, $num_c_semicolons)    = get_files_lines_semicolons_counts("c");
my ($num_f_files,   $num_f_lines, $num_f_semicolons)    = get_files_lines_semicolons_counts("f");
my ($num_F_files,   $num_F_lines, $num_F_semicolons)    = get_files_lines_semicolons_counts("F");

my $num_total_files = $num_hpp_files + $num_cpp_files + $num_H_files + $num_C_files + $num_h_files + $num_c_files + $num_f_files + $num_F_files;
my $num_total_lines = $num_hpp_lines + $num_cpp_lines + $num_H_lines + $num_C_lines + $num_h_lines + $num_c_lines + $num_f_lines + $num_F_lines;
my $num_total_semicolons = $num_hpp_semicolons + $num_cpp_semicolons + $num_H_semicolons + $num_C_semicolons + $num_h_semicolons
  + $num_c_semicolons + $num_f_semicolons + $num_F_semicolons;

my @pwd_list = split("/",$ENV{PWD});
print "\nStatistics for source files for ", $pwd_list[-1], " ...\n\n";

if($num_hpp_files) { print "Number of *.hpp files = $num_hpp_files\n"; }
if($num_cpp_files) { print "Number of *.cpp files = $num_cpp_files\n"; }
if($num_H_files)   { print "Number of *.H files   = $num_H_files\n"; }
if($num_C_files)   { print "Number of *.C files   = $num_C_files\n"; }
if($num_h_files)   { print "Number of *.h files   = $num_h_files\n"; }
if($num_c_files)   { print "Number of *.c files   = $num_c_files\n"; }
if($num_f_files)   { print "Number of *.f files   = $num_f_files\n"; }
if($num_F_files)   { print "Number of *.F files   = $num_F_files\n"; }
print "\n";
print "Total number of files = $num_total_files\n";

print "\n";

if($num_hpp_files) { print "Number of lines in *.hpp files = $num_hpp_lines\n"; }
if($num_cpp_files) { print "Number of lines in *.cpp files = $num_cpp_lines\n"; }
if($num_H_files)   { print "Number of lines in *.H files   = $num_H_lines\n"; }
if($num_C_files)   { print "Number of lines in *.C files   = $num_C_lines\n"; }
if($num_h_files)   { print "Number of lines in *.h files   = $num_h_lines\n"; }
if($num_c_files)   { print "Number of lines in *.c files   = $num_c_lines\n"; }
if($num_f_files)   { print "Number of lines in *.f files   = $num_f_lines\n"; }
if($num_F_files)   { print "Number of lines in *.F files   = $num_F_lines\n"; }
print "\n";
print "Total lines of code            = $num_total_lines\n";

print "\n";

if($num_hpp_files) { print "Number of semicolons in *.hpp files = $num_hpp_semicolons\n"; }
if($num_cpp_files) { print "Number of semicolons in *.cpp files = $num_cpp_semicolons\n"; }
if($num_H_files)   { print "Number of semicolons in *.H files   = $num_H_semicolons\n"; }
if($num_C_files)   { print "Number of semicolons in *.C files   = $num_C_semicolons\n"; }
if($num_h_files)   { print "Number of semicolons in *.h files   = $num_h_semicolons\n"; }
if($num_c_files)   { print "Number of semicolons in *.c files   = $num_c_semicolons\n"; }
#if($num_f_files)   { print "Number of semicolons in *.f files   = $num_f_semicolons\n"; }
#if($num_F_files)   { print "Number of semicolons in *.F files   = $num_F_semicolons\n"; }
print "\n";
print "Total number of semicolons          = $num_total_semicolons\n";

print "\n";

sub get_files_lines_semicolons_counts {
  my $file_ext = shift;
  my $output = `find $basedir -name "*.$file_ext" -exec wc -l '{}' ';'`;
  my @output_lines = split('\n',$output);
  my $num_files = 0;
  my $num_lines = 0;
  foreach(@output_lines) {
    ++$num_files;
    $num_lines += (split(" ",$_))[0];
#	  print "num_lines = $num_lines\n";
  }
  my $num_semicolons = 0;
  #$output = `find $basedir -name *.$file_ext -exec grep -Hc '\;' '{}' ';'`;
  $output = `find $basedir -name "*.$file_ext" -exec grep -Hc '\\;' '{}' ';'`;
  #print "output =\n$output\n";
  @output_lines = split('\n',$output);
  foreach(@output_lines) {
    $num_semicolons += (split(":",$_))[1];
#	  print "num_semicolons = $num_semicolons\n";
  }
  return ( $num_files, $num_lines, $num_semicolons );
}
