#!/bin/sh

# This script will replace all occurances of:
#
#    assert(anything)
#
# with
#
#    TEST_FOR_EXCEPT( !( anything ) ) 
#
# in the input file $1 in a faily safe way.
#
# This script is used as:
#
#    replace-assert-with-TEST_FOR_EXCEPT.sh filename
#
# and the file is modified in place.
#
# Note that you can run this script on a whole set of files using
# something like:
#
#    find . -name "*pp" -exec $ABS_DIR/replace-assert-with-TEST_FOR_EXCEPT.sh {} \;
#

sed -i "s/\([^a-zA-Z0-9_]\)assert(\(.\+\))/\1TEST_FOR_EXCEPT( \!( \2 ) )/g" $1
