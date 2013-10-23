#!/bin/sh

# run this from you Trilinos build directory 

# relative to Trilinos/packages/stk/stk_percept/build/build.dir

BASE_PATH=../../../../..
echo "BASE_PATH= " , `cd $BASE_PATH; pwd`

# Location of Trilinos source, where we build from

TRILINOS_CODE=$BASE_PATH

if [[ "$1" != "" ]] ; then
   TRILINOS_CODE=$1
fi
curwd=`pwd`

cd $TRILINOS_CODE/packages/PyTrilinos/doc/Doxygen/
make depend
make Percept_dox.i

cd $curwd
touch $TRILINOS_CODE/packages/PyTrilinos/src/stk/PerceptMesh.i
make -j8

pydoc -w ./packages/PyTrilinos/src/stk/PyPercept/PerceptMesh.py

cp PerceptMesh.html $TRILINOS_CODE/packages/stk/stk_percept/doc/



