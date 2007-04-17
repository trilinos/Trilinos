#!/bin/sh
#
# This script must be run from this directory
#

# Remove files from main so that they are not commited

rm ../configure*
rm ../Makefile.*
rm ../src/*.hpp
rm ../src/*.cpp
rm ../src/*.am
rm ../src/*.in
#rm ../doc/DoxyfileWeb
#rm ../doc/index.doc
#rm ../doc/images/Belos-Interfaces-Harder.gif

rm ../example/Makefile*
rm ../example/BlockGmres/Makefile*
rm ../example/BlockGmres/*.cpp
rm ../example/BlockGmres/*.hpp

rm ../test/Makefile*
rm ../test/BlockCG/Makefile*
rm ../test/BlockCG/*.cpp
rm ../test/BlockCG/*.hpp
rm ../test/BlockGmres/Makefile*
rm ../test/BlockGmres/*.cpp
rm ../test/BlockGmres/*.hpp

rm ../thyra/Makefile*.
rm ../thyra/src/Makefile* 
rm ../thyra/src/*.hpp 
rm ../thyra/example/Makefile*
rm ../thyra/example/LOWSFactory/Makefile*
rm ../thyra/example/LOWSFactory/Epetra/Makefile*
rm ../thyra/example/LOWSFactory/Epetra/*.cpp
rm ../thyra/example/LOWSFactory/Tpetra/Makefile*
rm ../thyra/example/LOWSFactory/Tpetra/*.cpp
rm ../thyra/example/LOWSFactory/Tpetra/*.hpp
rm ../thyra/test/Makefile*
rm ../thyra/test/LOWSFactory/Makefile*
rm ../thyra/test/LOWSFactory/*.hpp
rm ../thyra/test/LOWSFactory/*.cpp
rm ../thyra/test/MVOPTester/Makefile*
rm ../thyra/test/MVOPTester/*.cpp

