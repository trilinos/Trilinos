#!/bin/sh

cd teuchos/doc
./build.doxygen.documentation.sh
cd ../../
cd thyra/doc
./build.doxygen.documentation.sh
cd ../../
cd epetra/doc
./build.doxygen.documentation.sh
cd ../../
