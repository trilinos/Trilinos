#!/bin/bash

if [ $# -ne 1 ]; then
    echo "syntax: rebase.sh $TRILINOS_SRC/packages/muelu/test/interface"
    exit -1;
fi


if [ -d "default" ]; then
    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --linAlgebra=Tpetra
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor  --linAlgebra=Tpetra
    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --linAlgebra=Tpetra --heavytests
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --linAlgebra=Tpetra --heavytests
    ./MueLu_CreateOperator.exe --noKokkosRefactor --linAlgebra=Tpetra
    mpiexec -n 4 ./MueLu_CreateOperator.exe --noKokkosRefactor --linAlgebra=Tpetra

    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --linAlgebra=Epetra
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --linAlgebra=Epetra
    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --linAlgebra=Epetra --heavytests
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --linAlgebra=Epetra --heavytests
    ./MueLu_CreateOperator.exe --noKokkosRefactor --linAlgebra=Epetra
    mpiexec -n 4 ./MueLu_CreateOperator.exe --noKokkosRefactor --linAlgebra=Epetra

    pushd default/Output/
    source $1/default/Output/rebase.sh
    popd
else
    echo "Cannot rebase \"default\""
    exit 1
fi

if [ -d "kokkos" ]; then
    ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Tpetra --kokkos-num-threads=1
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Tpetra --kokkos-num-threads=1
    ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Tpetra --heavytests --kokkos-num-threads=1
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Tpetra --heavytests --kokkos-num-threads=1
    ./MueLu_CreateOperator.exe --useKokkosRefactor --linAlgebra=Tpetra --kokkos-num-threads=1
    mpiexec -n 4 ./MueLu_CreateOperator.exe --useKokkosRefactor --linAlgebra=Tpetra --kokkos-num-threads=1

    ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Epetra --kokkos-num-threads=1
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Epetra --kokkos-num-threads=1
    ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Epetra --heavytests --kokkos-num-threads=1
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --useKokkosRefactor --linAlgebra=Epetra --heavytests --kokkos-num-threads=1
    ./MueLu_CreateOperator.exe --useKokkosRefactor --linAlgebra=Epetra --kokkos-num-threads=1
    mpiexec -n 4 ./MueLu_CreateOperator.exe --useKokkosRefactor --linAlgebra=Epetra --kokkos-num-threads=1

    pushd kokkos/Output/
    source $1/kokkos/Output/rebase.sh
    popd
else
    echo "Cannot rebase \"kokkos\""
    exit 1
fi
