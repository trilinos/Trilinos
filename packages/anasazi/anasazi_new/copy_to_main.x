#!/bin/sh
#
# This script must be run from this directory
#

# Copy files from packages/anasazi/anasazi_new to packages/anasazi 

cp configure* ..
cp Makefile* ..

cp ./src/*.hpp ../src/.
cp ./src/*.cpp ../src/.
cp ./src/Makefile* ../src/.
cp ./src/*.in ../src/.

cp ./example/Makefile* ../example/.
mkdir ../example/BlockKrylovSchur
cp ./example/BlockKrylovSchur/*.cpp ../example/BlockKrylovSchur/.
cp ./example/BlockKrylovSchur/Makefile* ../example/BlockKrylovSchur/.

mkdir ../test
cp ./test/Makefile* ../test/.
mkdir ../test/ModalSolverUtils
cp ./test/ModalSolverUtils/*.cpp ../test/ModalSolverUtils/.
cp ./test/ModalSolverUtils/Makefile* ../test/ModalSolverUtils/.
mkdir ../test/BlockDavidson
cp ./test/BlockDavidson/*.cpp ../test/BlockDavidson/.
cp ./test/BlockDavidson/Makefile* ../test/BlockDavidson/.

mkdir ../util
cp ./util/Makefile* ../util/.
mkdir ../util/ModeLaplace
cp ./util/ModeLaplace/*.cpp ../util/ModeLaplace/.
cp ./util/ModeLaplace/*.h ../util/ModeLaplace/.
cp ./util/ModeLaplace/Makefile* ../util/ModeLaplace/.

