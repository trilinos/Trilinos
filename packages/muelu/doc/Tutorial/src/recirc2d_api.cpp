// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>

#include <Teuchos_StandardCatchMacros.hpp>

#include <Epetra_RowMatrix.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

// Galeri
#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>
#include <Galeri_Utils.h>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>

#include <MueLu_Utilities.hpp>

#include <MueLu_MutuallyExclusiveTime.hpp>

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

int main(int argc, char* argv[]) {
  return EXIT_SUCCESS;
}  // main
