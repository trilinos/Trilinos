#!/bin/sh
#
# This script must be run from this directory
#

# Remove files from main so that they are not commited

rm ../configure*
rm ../Makefile.*
rm ../src/*.hpp
rm ../src/*.cpp
rm ../src/*.am
rm ../src/*.in
rm ../doc/DoxyfileWeb
rm ../doc/index.doc
rm ../doc/images/Belos-Interfaces-Harder.gif

rm ../test/Makefile*
rm ../test/BlockCG/Makefile*
rm ../test/BlockCG/*.cpp
rm ../test/BlockCG/*.hpp
rm ../test/BlockGmres/Makefile*
rm ../test/BlockGmres/*.cpp
rm ../test/BlockGmres/*.hpp
