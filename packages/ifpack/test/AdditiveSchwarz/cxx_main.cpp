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
#include "Ifpack_AdditiveSchwarz.h"
#include "AztecOO.h"
#include "Ifpack_BlockPreconditioner.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#include "Ifpack_BlockJacobi.h"
#include "Ifpack_BlockGaussSeidel.h"
#include "Ifpack_Jacobi.h"
#include "Ifpack_GaussSeidel.h"
#include "Ifpack_SOR.h"
#include "Ifpack_SSOR.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_LinearPartitioner.h"
#include "Ifpack_GreedyPartitioner.h"
#include "Ifpack_METISPartitioner.h"

using namespace Trilinos_Util;

bool TestPreconditioner(Epetra_LinearProblem* Problem,
			int OverlapLevel = 0)
{

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());
  LHS.PutScalar(0.0);

  Teuchos::ParameterList List;
  Epetra_RowMatrix* A = Problem->GetMatrix();

  Ifpack_Preconditioner* Prec;
  
  Prec = new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(A,OverlapLevel);
  assert(Prec != 0);

  IFPACK_CHK_ERR(Prec->SetParameters(List));
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());

  // create the AztecOO solver
  AztecOO AztecOOSolver(*Problem);

  // specify solver
  AztecOOSolver.SetAztecOption(AZ_solver,AZ_gmres_condnum);
  AztecOOSolver.SetAztecOption(AZ_output,32);

  AztecOOSolver.SetPrecOperator(Prec);

  // solver. The solver should converge in one iteration,
  // or maximum two (numerical errors)
  AztecOOSolver.Iterate(1550,1e-8);

  double TrueResidual = AztecOOSolver.TrueResidual(); 
  // some output
  if( Problem->GetMatrix()->Comm().MyPID() == 0 ) {
    cout << "Solver performed " << AztecOOSolver.NumIters()
      << " iterations.\n";
    cout << "Norm of the true residual = " << TrueResidual << endl;
  }

  delete Prec;
  
  if (TrueResidual < 1e-5)
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
  const int NumPoints = 14400;

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // test the preconditioner
  int TestPassed = true;

  if (Comm.NumProc() == 1) {
    if (!TestPreconditioner(Problem,0))
      TestPassed = false;
  }
  else {
    for (int overlap = 0 ; overlap < 5 ; overlap++) {
      if (!TestPreconditioner(Problem,overlap))
	TestPassed = false;
    }
  }

  if (!TestPassed)
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif
  exit(EXIT_SUCCESS);
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
