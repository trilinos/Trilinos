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
rm ../example/BlockKrylovSchur/Makefile*
rm ../example/BlockKrylovSchur/*.cpp
rm -rf ../example/BlockKrylovSchur

rm ../test/Makefile*
rm ../test/ModalSolverUtils/Makefile*
rm ../test/ModalSolverUtils/*.cpp
rm -rf ../test/ModalSolverUtils
rm -rf ../test

rm ../util/Makefile*
rm ../util/ModeLaplace/Makefile*
rm ../util/ModeLaplace/*.cpp
rm ../util/ModeLaplace/*.h
rm -rf ../util/ModeLaplace
rm -rf ../util
