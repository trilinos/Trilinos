#!/bin/sh

# This script does the same thing as
# standard-C-CXX-header-upgrade-single-file.20070628.sh but does it
# recursively to all C++ source files in a directory tree.
#
# Warning! This script does not do replacements on *.h files since these could
# be C files.  If you need to replace on *.h files, please do this manually
# uisng the scirpt standard-C-CXX-header-upgrade-single-file.20070628.sh on a
# file by file basis.
#
# Warning!  Do not create a symbolic link to this script, just use a absolute
# or relative path from where it is defined in the Trilinos source tree.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

find . -name "*.hpp" -exec $_SCRIPT_DIR/standard-C-CXX-header-upgrade-single-file.20070628.sh {} \;
find . -name "*.cpp" -exec $_SCRIPT_DIR/standard-C-CXX-header-upgrade-single-file.20070628.sh {} \;
find . -name "*.H" -exec $_SCRIPT_DIR/standard-C-CXX-header-upgrade-single-file.20070628.sh {} \;
find . -name "*.C" -exec $_SCRIPT_DIR/standard-C-CXX-header-upgrade-single-file.20070628.sh {} \;
