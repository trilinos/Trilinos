#!/bin/sh
#
# This script must be run from this directory
#

# Remove files from main so that they are not commited

rm configure*

rm ../src/*.hpp
rm ../src/*.cpp
rm ../src/Makefile*

rm ../example/Makefile*
rm ../example/*.cpp
