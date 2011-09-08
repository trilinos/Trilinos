// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
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
    {
      cerr << "Caught exception: ";
      rhs.Print();
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
