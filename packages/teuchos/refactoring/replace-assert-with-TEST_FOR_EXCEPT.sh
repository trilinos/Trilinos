#!/bin/sh

# This script will replace all occurrences of:
#
#    assert(anything)
#
# with
#
#    TEUCHOS_TEST_FOR_EXCEPT( !( anything ) ) 
#
# in the input file $1 in a fairly safe way.
#
# This script is used as:
#
#    replace-assert-with-TEUCHOS_TEST_FOR_EXCEPT.sh filename
#
# and the file is modified in place.
#
# Note that you can run this script on a whole set of files using
# something like:
#
#    find . -name "*pp" -exec $ABS_DIR/replace-assert-with-TEUCHOS_TEST_FOR_EXCEPT.sh {} \;
#
#
# WARNING! This script does not work correctly for multi-line statements like:
#
#    assert(something); assert(anotherThing);
#
# as it will make this:
#
#    TEUCHOS_TEST_FOR_EXCEPT( !( something); assert(anotherThing ) );
#
# which is wrong wrong wrong.  Always, look over you code to make sure the right
# thing happened (or don't put multiple statements on the same line).
#
# WARNING! This script will not replace multi-line assert statements like:
#
#    assert(
#      something
#      );
#
# you must replace these manually!
#
# Also, don't forget to put in an:
#
#    #include "Teuchos_Assert.hpp"
#
# line in some header file that all of your code will include at least
# indirectly.
#

sed -i "s/\([^a-zA-Z0-9_]\)assert(\(.\+\))/\1TEUCHOS_TEST_FOR_EXCEPT( \!( \2 ) )/g" $1
