#!/bin/sh
#
# This script must be run from this directory
#

# Copy files from packages/belos/belos_new to packages/belos 

cp ./src/*.hpp ../src/.
cp ./src/*.cpp ../src/.
cp ./src/*.am ../src/.
cp ./src/*.in ../src/.

cp ./test/Makefile* ../test/.
cp ./test/BlockCG/Makefile* ../test/BlockCG/.
cp ./test/BlockCG/*.cpp ../test/BlockCG/.
cp ./test/BlockCG/*.hpp ../test/BlockCG/.
cp ./test/BlockGmres/Makefile* ../test/BlockGmres/.
cp ./test/BlockGmres/*.cpp ../test/BlockGmres/.
cp ./test/BlockGmres/*.hpp ../test/BlockGmres/.
