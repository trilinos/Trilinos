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
#include "Ifpack_METISPartitioner.h"

using namespace Trilinos_Util;

int TestPreconditioner(Epetra_LinearProblem* Problem,
		       int OverlapLevel, bool UseSparse)
{

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());
  LHS.PutScalar(0.0);

  Teuchos::ParameterList List;
  List.set("damping factor", .67);
  List.set("sweeps",5);
  List.set("local blocks", 4);
  List.set("overlap level", OverlapLevel);
  List.set("print level", 0);

  Epetra_RowMatrix* A = Problem->GetMatrix();
  Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);

  Ifpack_Preconditioner* Prec = 0;
  
  if (UseSparse)
    Prec = new Ifpack_CrsAdditiveSchwarz<Ifpack_BlockJacobi<Ifpack_SparseContainer<Ifpack_Amesos> > >(CrsA,OverlapLevel);
  else
    Prec = new Ifpack_CrsAdditiveSchwarz<Ifpack_BlockJacobi<Ifpack_DenseContainer> >(CrsA,OverlapLevel);

  assert(Prec != 0);

  IFPACK_CHK_ERR(Prec->SetParameters(List));
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
    return(0);
  else
    return(-1);

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
  const int NumPoints = 25;

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // test the preconditioner
  int TestPassed = true;
  for (int overlap = 1 ; overlap < 5 ; overlap++) {
    if (TestPreconditioner(Problem,overlap,true))
      TestPassed = false;
  }
  

#ifdef HAVE_MPI
  MPI_Finalize() ; 
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
