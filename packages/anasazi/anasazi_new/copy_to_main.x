#!/bin/sh
#
# This script must be run from this directory
#

# Copy files from packages/anasazi/anasazi_new to packages/anasazi 

cp configure* ..

cp ./src/*.hpp ../src/.
cp ./src/*.cpp ../src/.
cp ./src/Makefile* ../src/.

cp ./example/Makefile* ../example/.
mkdir ../example/ModeLaplace
cp ./example/ModeLaplace/*.cpp ../example/ModeLaplace/.
cp ./example/ModeLaplace/*.h ../example/ModeLaplace/.
cp ./example/ModeLaplace/Makefile* ../example/ModeLaplace/.
mkdir ../example/BlockKrylovSchur
cp ./example/BlockKrylovSchur/*.cpp ../example/BlockKrylovSchur/.
cp ./example/BlockKrylovSchur/Makefile* ../example/BlockKrylovSchur/.

