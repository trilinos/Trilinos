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

if( !(defined(@ARGV) && scalar(@ARGV)==1) ) {
  die "Error, this script takes one and only one argument for the base directory of a trilinos build.!\n";
}

# Takes just one options and that is the directory where Trilinos is built
my $TrilinosBuild = shift;

my $teuchos_makefile = "${TrilinosBuild}/packages/teuchos/src/Makefile";

my $cmnd =
"grep '^CC =' ${teuchos_makefile} ; "
."grep '^CXX =' ${teuchos_makefile} ; "
."grep '^libteuchos_a_AR =' ${teuchos_makefile} ; "
."grep '^CXXLD =' ${teuchos_makefile} ; "
."grep '^DEFS =' ${teuchos_makefile} ; "
."grep '^CPPFLAGS =' ${teuchos_makefile} ; "
."grep '^CFLAGS =' ${teuchos_makefile} ; "
."grep '^CXXFLAGS =' ${teuchos_makefile} ; "
."grep '^FFLAGS =' ${teuchos_makefile} ; "
."grep '^LDFLAGS =' ${teuchos_makefile} ; "
."grep '^FLIBS =' ${teuchos_makefile} ; "
."grep '^LIBS =' ${teuchos_makefile} ; "
."grep '^BLAS_LIBS =' ${teuchos_makefile} ; "
."grep '^LAPACK_LIBS =' ${teuchos_makefile} ; "
."grep '^prefix =' ${teuchos_makefile} ; "
."grep '^RANLIB =' ${teuchos_makefile} "
;

#print $cmnd;

my $sys_return = system($cmnd);

#print "\nsys_return = ${sys_return}\n";

exit($sys_return==0 ? 0 : 1);
