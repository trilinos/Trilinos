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

my $success = 1;
my $result;
$result = system ('./RefCountPtr/RefCountPtr_test.exe');
$success = 0 if (!$result);
exit ($success ? 0 : -1 );
