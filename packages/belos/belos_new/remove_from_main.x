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
rm ../src/interfaces/*.hpp
rmdir ../src/interfaces
rm ../src/implementations/*.hpp
rmdir ../src/implementations
rm ../doc/DoxyfileWeb
rm ../doc/index.doc

rm ../test/Makefile*
rm ../test/BlockCG/Makefile*
rm ../test/BlockCG/*.cpp
rm ../test/BlockCG/*.hpp
rm ../test/BlockGmres/Makefile*
rm ../test/BlockGmres/*.cpp
rm ../test/BlockGmres/*.hpp
rm ../test/NonblockGmres/Makefile*
rm ../test/NonblockGmres/*.cpp
rm ../test/NonblockGmres/*.hpp
rm ../test/NonblockGmres/*.pl
rmdir ../test/NonblockGmres
rm ../test/NonblockCg/Makefile*
rm ../test/NonblockCg/*.cpp
#rm ../test/NonblockCg/*.hpp
rm ../test/NonblockCg/*.pl
rmdir ../test/NonblockCg
