#!/usr/bin/perl -w
#
# This perl script graps a bunch of make macro definitions
# generated for Teuchos that can be used in other makefiles.
# This is dumped to stdout and can be redirected to build
# a makefile.
#
# Note, this script must be maintained to be current for
# the Teuchos makefile.
#
use strict;

if( !(defined(@ARGV) && scalar(@ARGV)==2) ) {
  die "Error, this script takes two and only two arguments (makefile_name package_name).!\n";
}

my $makefile_name = shift;
my $package_name  = shift;

#
# List the macros you want to grep and include in the output
#
my @macros =
	(
	 "CC"
	 ,"CXX"
	 ,"F77"
	 ,"CXXLD"
	 ,"DEFS"
	 ,"CPPFLAGS"
	 ,"CFLAGS"
	 ,"CXXFLAGS"
	 ,"FFLAGS"
	 ,"LDFLAGS"
	 ,"FLIBS"
	 ,"BLAS_LIBS"
	 ,"LAPACK_LIBS"
	 ,"prefix"
	 ,"AR"
	 ,"ALTERNATE_AR"
	 ,"libteuchos_a_AR"
	 ,"RANLIB"
	 );

open FILE_IN, "<$makefile_name" || die "The file $makefile_name could not be opended for input\n";
my @makefile_name_array = <FILE_IN>;
close FILE_IN;

#
# Find the above macros and append "${package_name}_" to the beginning.
#
my @new_macros;
my $add_next_line = 0;
foreach( @makefile_name_array ) {
	my $line = $_;
	if($add_next_line) {
		push @new_macros, $line;
		if( substr($line,-1,1) eq "\\" ) {
			$add_next_line = 1;
		}
		else {
			$add_next_line = 0;
		}
		next;
	}
	#print "Line = $line";
	foreach( @macros ) {
		my $macro_search = "^${_} ";
		#print "Macro search = \'$macro_search\'\n";
		if( $line=~/$macro_search/ ) {
			#print "Adding Macro!\n";
      my $find_str = '\(CXX\)';
      my $replace_str = "(${package_name}_CXX)";
      $line=~s/$find_str/$replace_str/;
			push @new_macros, "${package_name}_${line}";
			if( substr($line,-2,1) eq "\\" ) {
				$add_next_line = 1;
			}
			else {
				$add_next_line = 0;
			}
		}
	}
}

print join("",@new_macros);
