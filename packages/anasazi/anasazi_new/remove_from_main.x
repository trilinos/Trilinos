#!/bin/sh
#
# This script must be run from this directory
#

# Remove files from main so that they are not commited

rm ../configure*

rm ../src/*.hpp
rm ../src/*.cpp
rm ../src/Makefile*

rm ../example/Makefile*
rm ../example/ModeLaplace/Makefile*
rm ../example/ModeLaplace/*.cpp
rm ../example/ModeLaplace/*.h
rm -rf ../example/ModeLaplace
rm ../example/BlockArnoldi/Makefile*
rm ../example/BlockArnoldi/*.cpp
rm -rf ../example/BlockArnoldi
