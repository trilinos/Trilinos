#!/bin/sh
#
# This script must be run from this directory
#

# Copy files from packages/belos to packages/anasazi/anasazi_new

cp ../configure* .

cp ../src/*.hpp ./src/.
cp ../src/*.cpp ./src/.
cp ../src/Makefile* ./src/.

cp ../example/Makefile* ./example/.
cp ../example/*.cpp ./example/.

# Remove files from main so that they are not commited
./remove_from_main.x
