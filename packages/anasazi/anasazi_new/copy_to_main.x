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
mkdir ../example/ModeLaplace
cp ./example/ModeLaplace/*.cpp ../example/ModeLaplace/.
cp ./example/ModeLaplace/*.h ../example/ModeLaplace/.
cp ./example/ModeLaplace/Makefile* ../example/ModeLaplace/.
mkdir ../example/BlockKrylovSchur
cp ./example/BlockKrylovSchur/*.cpp ../example/BlockKrylovSchur/.
cp ./example/BlockKrylovSchur/Makefile* ../example/BlockKrylovSchur/.

mkdir ../test
cp ./test/Makefile* ../test/.
mkdir ../test/ModalSolverUtils
cp ./test/ModalSolverUtils/*.cpp ../test/ModalSolverUtils/.
cp ./test/ModalSolverUtils/Makefile* ../test/ModalSolverUtils/.

