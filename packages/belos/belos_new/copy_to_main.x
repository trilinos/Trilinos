#!/bin/sh
#
# This script must be run from this directory
#

# Copy files from packages/belos/belos_new to packages/belos 

cp ./configure* ../.
cp ./Makefile.* ../.
cp ./src/*.h ../src/.
cp ./src/*.hpp ../src/.
cp ./src/*.cpp ../src/.
cp ./src/*.am ../src/.
cp ./src/*.in ../src/.

#cp ./doc/DoxyfileWeb ../doc/.
#cp ./doc/index.doc ../doc/.
#cp ./doc/images/Belos-Interfaces-Harder.gif ../doc/images/.

cp ./example/Makefile* ../example/.
cp ./example/BlockGmres/Makefile* ../example/BlockGmres/.
cp ./example/BlockGmres/*.cpp ../example/BlockGmres/.
cp ./example/BlockCG/Makefile* ../example/BlockCG/.
cp ./example/BlockCG/*.cpp ../example/BlockCG/.

cp ./test/Makefile* ../test/.
cp ./test/BlockCG/Makefile* ../test/BlockCG/.
cp ./test/BlockCG/*.cpp ../test/BlockCG/.
cp ./test/BlockCG/*.hpp ../test/BlockCG/.
cp ./test/BlockGmres/Makefile* ../test/BlockGmres/.
cp ./test/BlockGmres/*.cpp ../test/BlockGmres/.
cp ./test/BlockGmres/*.hpp ../test/BlockGmres/.

cp ./thyra/Makefile* ../thyra/.
cp ./thyra/src/Makefile* ../thyra/src/.
cp ./thyra/src/*.hpp ../thyra/src/.
cp ./thyra/example/Makefile* ../thyra/example/.
cp ./thyra/example/LOWSFactory/Makefile* ../thyra/example/LOWSFactory/.
cp ./thyra/example/LOWSFactory/Epetra/Makefile* ../thyra/example/LOWSFactory/Epetra/.
cp ./thyra/example/LOWSFactory/Epetra/*.cpp ../thyra/example/LOWSFactory/Epetra/.
cp ./thyra/example/LOWSFactory/Tpetra/Makefile* ../thyra/example/LOWSFactory/Tpetra/.
cp ./thyra/example/LOWSFactory/Tpetra/*.cpp ../thyra/example/LOWSFactory/Tpetra/.
cp ./thyra/example/LOWSFactory/Tpetra/*.hpp ../thyra/example/LOWSFactory/Tpetra/.
cp ./thyra/test/Makefile* ../thyra/test/.
cp ./thyra/test/LOWSFactory/Makefile* ../thyra/test/LOWSFactory/.
cp ./thyra/test/LOWSFactory/*.hpp ../thyra/test/LOWSFactory/.
cp ./thyra/test/LOWSFactory/*.cpp ../thyra/test/LOWSFactory/.
cp ./thyra/test/MVOPTester/Makefile* ../thyra/test/MVOPTester/.
cp ./thyra/test/MVOPTester/*.cpp ../thyra/test/MVOPTester/.

