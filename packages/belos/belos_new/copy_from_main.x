#!/bin/sh
#
# This script must be run from this directory
#

# Copy files from packages/belos to packages/belos/belos_new

cp ../src/*.hpp ./src/.
cp ../src/*.cpp ./src/.
cp ../src/*.am ./src/.
cp ../src/*.in ./src/.

cp ../test/Makefile* ./test/.
cp ../test/BlockCG/Makefile* ./test/BlockCG/.
cp ../test/BlockCG/*.cpp ./test/BlockCG/.
cp ../test/BlockGmres/Makefile* ./test/BlockGmres/.
cp ../test/BlockGmres/*.cpp ./test/BlockGmres/.

# Remove files from main so that they are not commited

rm ../src/*.hpp
rm ../src/*.cpp
rm ../src/*.am
rm ../src/*.in

rm ../test/Makefile*
rm ../test/BlockCG/Makefile*
rm ../test/BlockCG/*.cpp
rm ../test/BlockGmres/Makefile*
rm ../test/BlockGmres/*.cpp
