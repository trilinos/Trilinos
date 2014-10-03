#!/bin/bash

_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*remove_std_tribits_includes_r.sh/\1/g"`
#echo $_SCRIPT_DIR

echo
echo "Replacing all standard TriBITS includes in CMakeLists.txt files ..."
echo

find . -name CMakeLists.txt -exec $_SCRIPT_DIR/remove_std_tribits_includes.py {} {} \;

echo
echo "Replacing all standard TriBITS includes in CMakeLists.tribits files ..."
echo

find . -name CMakeLists.tribits -exec $_SCRIPT_DIR/remove_std_tribits_includes.py {} {} \;

echo
echo "Replacing all standard TriBITS includes in *.cmake files ..."
echo

find . -name "*.cmake" -exec $_SCRIPT_DIR/remove_std_tribits_includes.py {} {} \;
