// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
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

int main(int argv, char* argc[])
{
#ifdef HAVE_MPI
  MPI_Init(&argv, &argc);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Here we create the linear problem
  //
  //   Matrix * LHS = RHS
  //
  // with Matrix arising from a 5-point formula discretization.
  
  Epetra_Map*         Map = 0;
  Epetra_RowMatrix*   Matrix = 0;

  Teuchos::ParameterList GaleriList;
  // dimension of the problem is nx x ny
  GaleriList.set("nx", 10 * Comm.NumProc());
  GaleriList.set("ny", 10);
  // total number of processors is mx x my
  GaleriList.set("mx", Comm.NumProc());
  GaleriList.set("my", 1);

  try
  {
    Map = CreateMap("Cartesian2D", Comm, GaleriList);
    Matrix = CreateCrsMatrix("Laplace2D", Map, GaleriList);
    Epetra_Vector ExactSolution(*Map); ExactSolution.Random();
    Epetra_Vector LHS(*Map); LHS.PutScalar(0.0);
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
  }
  catch (Galeri::Exception& rhs)
  {
    if (Comm.MyPID() == 0)
      rhs.Print();
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
