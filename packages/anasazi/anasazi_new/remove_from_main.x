#!/bin/sh
#
# This script must be run from this directory
#

# Remove files from main so that they are not commited

rm ../configure*
rm ../Makefile*

rm ../src/*.hpp
rm ../src/*.cpp
rm ../src/Makefile*
rm ../src/*.in

rm ../example/Makefile*
rm ../example/ModeLaplace/Makefile*
rm ../example/ModeLaplace/*.cpp
rm ../example/ModeLaplace/*.h
rm -rf ../example/ModeLaplace
rm ../example/BlockKrylovSchur/Makefile*
rm ../example/BlockKrylovSchur/*.cpp
rm -rf ../example/BlockKrylovSchur

rm ../test/Makefile*
rm ../test/ModalSolverUtils/Makefile*
rm ../test/ModalSolverUtils/*.cpp
rm -rf ../test/ModalSolverUtils
rm -rf ../test
