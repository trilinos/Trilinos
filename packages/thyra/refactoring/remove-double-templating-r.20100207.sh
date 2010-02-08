#!/usr/bin/env sh
#
# This script removes the last remnants of double templating
# <RangeScalar,DomainScalar> and replaces it with just <Scalar>.
# 
# Run this script from a base of directories.
#

_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/remove-double-templating.*/\1/g"`
echo $_SCRIPT_DIR

find . -name "*pp" -exec $_SCRIPT_DIR/remove-double-templating.sh {} \;
