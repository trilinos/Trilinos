#!/bin/sh
#
# This script must be run from this directory
#

# Copy files from packages/belos/belos_new to packages/belos 

cp ./configure* ../.
cp ./Makefile.* ../.
cp ./src/*.hpp ../src/.
cp ./src/*.cpp ../src/.
cp ./src/*.am ../src/.
cp ./src/*.in ../src/.
mkdir ../src/interfaces
cp ./src/interfaces/*.hpp ../src/interfaces/.
mkdir ../src/implementations
cp ./src/implementations/*.hpp ../src/implementations/.
cp ./doc/DoxyfileWeb ../doc/.
cp ./doc/index.doc ../doc/.
cp ./doc/images/Belos-Interfaces-Harder.gif ../doc/images/.

cp ./test/Makefile* ../test/.
cp ./test/BlockCG/Makefile* ../test/BlockCG/.
cp ./test/BlockCG/*.cpp ../test/BlockCG/.
cp ./test/BlockCG/*.hpp ../test/BlockCG/.
cp ./test/BlockGmres/Makefile* ../test/BlockGmres/.
cp ./test/BlockGmres/*.cpp ../test/BlockGmres/.
cp ./test/BlockGmres/*.hpp ../test/BlockGmres/.
mkdir ../test/NonblockGmres
cp ./test/NonblockGmres/Makefile* ../test/NonblockGmres/.
cp ./test/NonblockGmres/*.cpp ../test/NonblockGmres/.
cp ./test/NonblockGmres/*.hpp ../test/NonblockGmres/.
cp ./test/NonblockGmres/*.pl ../test/NonblockGmres/.
mkdir ../test/NonblockCg
cp ./test/NonblockCg/Makefile* ../test/NonblockCg/.
cp ./test/NonblockCg/*.cpp ../test/NonblockCg/.
#cp ./test/NonblockCg/*.hpp ../test/NonblockCg/.
cp ./test/NonblockCg/*.pl ../test/NonblockCg/.
