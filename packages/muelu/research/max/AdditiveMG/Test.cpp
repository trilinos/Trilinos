// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"

// MueLu
#include "MueLu.hpp"
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

// Belos
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"

using namespace Galeri;

// =========== //
// main driver //
// =========== //

int main(int argc, char* argv[]) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Here we create the linear problem
  //
  //   Matrix * LHS = RHS
  //
  // with Matrix arising from a 5-point formula discretization.

  Epetra_Map* Map          = 0;
  Epetra_RowMatrix* Matrix = 0;

  Teuchos::ParameterList GaleriList;
  // dimension of the problem is nx x ny
  GaleriList.set("nx", 10 * Comm.NumProc());
  GaleriList.set("ny", 10);
  // total number of processors is mx x my
  GaleriList.set("mx", Comm.NumProc());
  GaleriList.set("my", 1);

  try {
    Map    = CreateMap("Cartesian2D", Comm, GaleriList);
    Matrix = CreateCrsMatrix("Laplace2D", Map, GaleriList);
    Epetra_Vector ExactSolution(*Map);
    ExactSolution.Random();
    Epetra_Vector LHS(*Map);
    LHS.PutScalar(0.0);
    Epetra_Vector RHS(*Map);

    Matrix->Multiply(false, ExactSolution, RHS);

    Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

    // at this point any object that understand Epetra_LinearProblem can be
    // used, for example AztecOO, Amesos. IFPACK and ML can be used to define a
    // preconditioner for Matrix. Here we use a simple solver, based on
    // LAPACK, that is meant for simple testing only.

    Solve(Problem);

    // and we compute the norm of the true residual.
    double ResidualNorm = ComputeNorm(Matrix, &LHS, &RHS);

    if (Comm.MyPID() == 0)
      cout << ResidualNorm << endl;

    delete Map;
    delete Matrix;
  } catch (Galeri::Exception& rhs) {
    if (Comm.MyPID() == 0)
      rhs.Print();
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}
