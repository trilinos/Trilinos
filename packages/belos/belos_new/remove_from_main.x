#!/bin/sh
#
# This script must be run from this directory
#

# Remove files from main so that they are not commited

rm ../src/*.hpp
rm ../src/*.cpp
rm ../src/*.am
rm ../src/*.in

rm ../test/Makefile*
rm ../test/BlockCG/Makefile*
rm ../test/BlockCG/*.cpp
rm ../test/BlockCG/*.hpp
rm ../test/BlockGmres/Makefile*
rm ../test/BlockGmres/*.cpp
rm ../test/BlockGmres/*.hpp
