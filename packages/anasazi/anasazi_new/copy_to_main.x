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
cp ./example/*.cpp ../example/.
