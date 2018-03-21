#!/bin/env bash

#The purpose of this file is to generate a list of package enables that
#when used with forward packages turned on will result in the right set of
#tests being run to test a pull request. The current implementation is neither
#maintainable or precise. It is a proof-of-concept at this point.

#file with list of changed files - generated currently with 
#git diff origin/$TRILINOS_TARGET_BRANCH --numstat > ../gitchanges.txt
#assumed to be in the directory where this script is executed from

declare PackageEnables=""

#This could be improved. There are too many cases where you actually
#don't need to build anything if something changes in the cmake directory.
#Right now the goal is to build everything for changes to general cmake
#files, such as those in Trilinos/cmake/, and Trilinos/CTestConfig.cmake,
#Trilinos/PackagesList.cmake, etc. as changes to those files can have
#impact across packages. However, if files change in
#Trilinos/packages/<package>/cmake, we will just build <package> and its
#forward dependencies.

if grep -i cmake gitchanges.txt |grep -qv packages; then
      PackageEnables+="-DTrilinos_ENABLE_ALL_PACKAGES=ON "
else #If we aren't building everything, figure out which packages to bulid

  if grep -q commonTools/gtest/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Gtest=ON "
  fi

  if grep -q packages/ThreadPool/ gitchanges.txt; then
	PackageEnables+="-DTrilinos_ENABLE_ThreadPool=ON "
  fi

  if grep -q packages/kokkos/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Kokkos=ON "
  fi

  if grep -q packages/teuchos/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Teuchos=ON "
  fi


  if grep -q packages/kokkos-kernels/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_KokkosKernels=ON "
  fi

  if grep -q packages/rtop/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_RTOp=ON "
  fi

  if grep -q packages/sacado/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Sacado=ON "
  fi

  if grep -q packages/minitensor/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_MiniTensor=ON "
  fi

  if grep -q packages/epetra/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Epetra=ON "
  fi

  if grep -q packages/zoltan/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Zoltan=ON "
  fi

  if grep -q packages/shards/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Shards=ON "
  fi

  if grep -q packages/globipack/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_GlobiPack=ON "
  fi

  if grep -q packages/triutils/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Triutils=ON "
  fi

  if grep -q packages/tpetra/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Tpetra=ON "
  fi

  if grep -q packages/common/auxiliarySoftware/SuiteSparse/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_TrilinosSS=ON "
  fi

  if grep -q packages/epetraext/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_EpetraExt=ON "
  fi

  if grep -q packages/domi/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Domi=ON "
  fi

  if grep -q packages/thyra/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Thyra=ON "
  fi

  if grep -q packages/xpetra/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Xpetra=ON "
  fi

  if grep -q packages/optipack/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_OptiPack=ON "
  fi

  if grep -q packages/isorropia/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Isorropia=ON "
  fi

  if grep -q packages/aztecoo/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_AztecOO=ON "
  fi

  if grep -q packages/galeri/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Galeri=ON "
  fi

  if grep -q packages/amesos/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Amesos=ON "
  fi

  if grep -q packages/pamgen/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Pamgen=ON "
  fi

  if grep -q packages/zoltan2/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Zoltan2=ON "
  fi

  if grep -q packages/ifpack/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Ifpack=ON "
 fi

  if grep -q packages/ml/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_ML=ON "
  fi

  if grep -q packages/belos/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Belos=ON "
  fi

  if grep -q packages/shylu/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_ShyLU=ON "
  fi

  if grep -q packages/amesos2/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Amesos2=ON "
  fi

  if grep -q packages/seacas/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_SEACAS=ON "
  fi

  if grep -q packages/anasazi/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Anasazi=ON "
  fi

  if grep -q packages/ifpack2/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Ifpack2=ON "
  fi

  if grep -q packages/stratimikos/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Stratimikos=ON "
  fi

  if grep -q packages/fei/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_FEI=ON "
  fi

  if grep -q packages/teko/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Teko=ON "
  fi

  if grep -q packages/intrepid/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Intrepid=ON "
  fi

  if grep -q packages/intrepid2/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Intrepid2=ON "
  fi

  if grep -q packages/stk/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_STK=ON "
  fi

  if grep -q packages/phalanx/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Phalanx=ON "
  fi

  if grep -q packages/nox/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_NOX=ON "
  fi

  if grep -q packages/muelu/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_MueLu=ON "
  fi

  if grep -q packages/rythmos/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Rythmos=ON "
  fi

  if grep -q packages/tempus/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Tempus=ON "
  fi

  if grep -q packages/stokhos/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Stokhos=ON "
  fi

  if grep -q packages/rol/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_ROL=ON "
  fi

  if grep -q packages/piro/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Piro=ON "
  fi

  if grep -q packages/panzer/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Panzer=ON "
  fi

  if grep -q packages/trilinoscouplings/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_TrilinosCouplings=ON "
  fi

  if grep -q packages/pike/ gitchanges.txt; then
        PackageEnables+="-DTrilinos_ENABLE_Pike=ON "
  fi
fi

#Turn some things off that are currently problematic
#PackageEnables+="-DIfpack2_ENABLE_TESTS=OFF -DMueLu_ENABLE_TESTS=OFF -DAnasazi_ENABLE_TESTS=OFF "

echo "$PackageEnables"

