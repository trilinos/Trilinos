#!/usr/bin/perl -w
#
# Must be run from the base directory

use FileHandle;
use Cwd;

my $basedir = Cwd::cwd();

my ($num_hpp_files,  $num_hpp_lines) = get_files_lines_counts("hpp");
my ($num_cpp_files, $num_cpp_lines)  = get_files_lines_counts("cpp");
my ($num_H_files,   $num_H_lines)    = get_files_lines_counts("H");
my ($num_C_files,   $num_C_lines)    = get_files_lines_counts("C");
my ($num_h_files,   $num_h_lines)    = get_files_lines_counts("h");
my ($num_c_files,   $num_c_lines)    = get_files_lines_counts("c");
my ($num_f_files,   $num_f_lines)    = get_files_lines_counts("f");
my ($num_F_files,   $num_F_lines)    = get_files_lines_counts("F");

my $num_total_files = $num_hpp_files + $num_cpp_files + $num_H_files + $num_C_files + $num_h_files + $num_c_files + $num_f_files + $num_F_files;
my $num_total_lines = $num_hpp_lines + $num_cpp_lines + $num_H_lines + $num_C_lines + $num_h_lines + $num_c_lines + $num_f_lines + $num_F_lines;

print "\nStatistics for source files ...\n\n";

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
if($num_H_files)   { print "Number of lines in *.H files   = $num_hpp_lines\n"; }
if($num_C_files)   { print "Number of lines in *.C files   = $num_cpp_lines\n"; }
if($num_h_files)   { print "Number of lines in *.h files   = $num_h_lines\n"; }
if($num_c_files)   { print "Number of lines in *.c files   = $num_c_lines\n"; }
if($num_f_files)   { print "Number of lines in *.f files   = $num_f_lines\n"; }
if($num_F_files)   { print "Number of lines in *.F files   = $num_F_lines\n"; }
print "\n";
print "Total lines of code            = $num_total_lines\n";

print "\n";

sub get_files_lines_counts {
  my $file_ext = shift;
  my $output = `find $basedir -name *.$file_ext -exec wc -l '{}' ';'`;
  return count_lines_files_from_output(split('\n',$output));
}

sub count_lines_files_from_output {
  my @output_lines = @_;
  my $num_files = 0;
  my $num_lines = 0;
  foreach(@output_lines) {
	++$num_files;
	$num_lines += (split(" ",$_))[0];
#	print "num_lines = $num_lines\n";
  }
  return ( $num_files, $num_lines );
}
