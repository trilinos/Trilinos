#!/usr/bin/perl -w
use strict;
use strict 'refs';
#
# RAB: 2003/11/14: According the the documentation that I can find
# it looks like this script should be invoked from:
#
#    BUILD_DIR/packages/teuchos/test
#
# and therefore I will write all of my relative paths from that
#
printf
  "\n******************************".
  "\n*** Running Teuchos tests ****".
  "\n******************************\n";

my $success = 1;  # Boolean (false=0,true=nonzero)
my $result;       # success=0, failure=nonzero
$result = system ('./RefCountPtr/RefCountPtr_test.exe --quiet');
$success = 0 if ($result != 0);
$result = system ('./BLAS/BLAS_test.exe');
$success = 0 if ($result != 0);
exit ($success ? 0 : -1 );
