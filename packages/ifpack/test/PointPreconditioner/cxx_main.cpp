// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AZTECOO) && defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_Jacobi.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_CrsAdditiveSchwarz.h"
#include "Ifpack_BlockGaussSeidel.h"
#include "Ifpack_Jacobi.h"
#include "Ifpack_GaussSeidel.h"
#include "Ifpack_SOR.h"
#include "Ifpack_SSOR.h"

using namespace Trilinos_Util;

bool TestPreconditioner(string PrecType,
		       CrsMatrixGallery& Gallery)
{

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());

  // Set up the list
  Teuchos::ParameterList List;
  List.set("point: damping factor", 1.0);
  List.set("point: sweeps",1550);
  List.set("point: print frequency", 1000);

  Ifpack_Preconditioner* PointPrec = 0;

  if (PrecType == "Jacobi")
    PointPrec = new Ifpack_Jacobi(A);
  else if (PrecType == "Gauss-Seidel")
    PointPrec = new Ifpack_GaussSeidel(A);
  else if (PrecType == "SOR")
    PointPrec = new Ifpack_SOR(A);
  else if (PrecType == "SSOR")
    PointPrec = new Ifpack_SSOR(A);

  assert (PointPrec != 0);

  PointPrec->SetParameters(List);
  PointPrec->Compute();
  // use the preconditioner as solver, with 1550 iterations
  PointPrec->ApplyInverse(RHS,LHS);

  delete PointPrec;

  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);

  if( A->Comm().MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }
  
  if (residual < 1e-5) 
    return(true);
  else
    return(false);
}

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // size of the global matrix. 
  const int NumPoints = 900;

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");


  // test the preconditioner
  int TestPassed = true;

  if (TestPreconditioner("Jacobi",Gallery))
    TestPassed = false;

  if (TestPreconditioner("Gauss-Seidel",Gallery))
    TestPassed = false;

  if (TestPreconditioner("SOR",Gallery))
    TestPassed = false;

  if (TestPreconditioner("SSOR",Gallery))
    TestPassed = false;

#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  if (TestPassed)
    exit(EXIT_SUCCESS);
  else
    exit(EXIT_FAILURE);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  puts("please configure IFPACK with --eanble-aztecoo --enable-teuchos");
  puts("--enable-amesos to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
