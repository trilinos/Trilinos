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
#include "Ifpack_Preconditioner.h"
#include "Ifpack.h"
#include "AztecOO.h"
#include "Ifpack_CrsIct.h"
#include "Ifpack_DropFilter.h"
#include "Amesos_TestRowMatrix.h"

using namespace Trilinos_Util;

int TestPreconditioner(string PrecType, Epetra_LinearProblem* Problem,
		       int OverlapLevel = 0)
{

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());
  LHS.PutScalar(0.0);

  Teuchos::ParameterList List;
  Epetra_RowMatrix* A = Problem->GetMatrix();
  Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);

  Amesos_TestRowMatrix RowMatrixA(A);

#define NEW_FACT
#ifdef NEW_FACT
  Ifpack Factory;

  int MaxRowEntries = 3;
  double DropValue = 0.0;
  double MinAbsDiagValue = 0.0;
  double AddToDiagValue = 0.0;
  /*
  Ifpack_DropFilter Adrop(&RowMatrixA, DropValue, MinAbsDiagValue,
			  AddToDiagValue, MaxRowEntries);
			  */
  Ifpack_Preconditioner* Prec = Factory.Create(PrecType, A, OverlapLevel);

  assert(Prec != 0);

  List.set("fact: level-of-fill", OverlapLevel);
  List.set("partitioner: local parts", 74);
  List.set("partitioner: type", "metis");
  List.set("partitioner: overlap", 0);
  List.set("schwarz: use RCM reordering", true);
  IFPACK_CHK_ERR(Prec->SetParameters(List));
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
#else
   int    LevelFill = OverlapLevel;
   double DropTol = 0.0;
   double Condest;

      Ifpack_CrsIct * Prec = NULL;
       Prec = new Ifpack_CrsIct(*CrsA,DropTol,LevelFill);
          Prec->InitValues(*CrsA);
       Prec->Factor();
      Prec->Condest(false,Condest);
               cout << Condest << endl;
#endif

  // create the AztecOO solver
  AztecOO AztecOOSolver(*Problem);

  // specify solver
  AztecOOSolver.SetAztecOption(AZ_solver,AZ_gmres);
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
  
  exit(0);
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
  const int NumPoints = 10000;

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // test the preconditioner
  int TestPassed = true;
  for (int overlap = 0 ; overlap < 5 ; overlap++) {
    if (TestPreconditioner("block relaxation",Problem,overlap))
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
