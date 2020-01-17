#!/bin/bash

if [ $# -ne 1 ]; then
    echo "syntax: rebase.sh $TRILINOS_SRC/packages/muelu/test/interface"
    exit -1;
fi


if [ -d "default" ]; then
    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra
    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra --heavytests
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra --heavytests
    ./MueLu_CreateOperator.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra
    mpiexec -n 4 ./MueLu_CreateOperator.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra

    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra
    ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra --heavytests
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra --heavytests
    ./MueLu_CreateOperator.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra
    mpiexec -n 4 ./MueLu_CreateOperator.exe --noKokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra

    pushd default/Output/
    source $1/default/Output/rebase.sh
    popd
else
    echo "Cannot rebase \"default\""
    exit 1
fi

if [ -d "kokkos" ]; then
    ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra
    ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra --heavytests
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra --heavytests
    ./MueLu_CreateOperator.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra
    mpiexec -n 4 ./MueLu_CreateOperator.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Tpetra

    ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra
    ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra --heavytests
    mpiexec -n 4 ./MueLu_ParameterListInterpreter.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra --heavytests
    ./MueLu_CreateOperator.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra
    mpiexec -n 4 ./MueLu_CreateOperator.exe --kokkosRefactor --instantiation="DOUBLE_INT_INT" --linAlgebra=Epetra

    pushd kokkos/Output/
    source $1/kokkos/Output/rebase.sh
    popd
else
    echo "Cannot rebase \"kokkos\""
    exit 1
fi
