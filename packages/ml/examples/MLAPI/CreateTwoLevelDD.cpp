
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_include.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_TwoLevelDDHybrid.h"
#include "MLAPI_TwoLevelDDAdditive.h"
#include "MLAPI_EpetraPreconditioner.h"
#include "MLAPI_Smoother.h"

using namespace Teuchos;
using namespace Trilinos_Util;
using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  // initialize the MLAPI workspace
  Init();
#endif

  try {

    int XNodes = 10;
    int YNodes = 10;
    // define a very simple linear map
    int NumGlobalElements = XNodes * YNodes;
    int NumMyElements = NumGlobalElements / NumProc();
    if (MyPID() == NumProc() - 1)
      NumMyElements += NumGlobalElements % NumProc();

    Teuchos::ParameterList MLList;

    Space FineSpace(NumMyElements);
    Operator FineMatrix;

    cerr << "TO BE FIXED..." << endl;
    exit(0);

    // define the aggregation-based prolongator operator. This
    // is the non-smoothed prolongator (tentative prolongator)
    SetPrintLevel(10);
    MLList.set("aggregation: type","Uncoupled");
    MLList.set("aggregation: nodes per aggregate",4);
    Operator Ptent = BuildP(FineMatrix,MLList);
    // get the diagonal of the fine matrix
    DoubleVector Diag = Diagonal(FineMatrix);
    // invert the elements on Diag
    Diag.Reciprocal();
    // build a matrix which contains only diagonal elements
    // given by vector Diag
    Operator Dinv = Diagonal(FineSpace,FineSpace,Diag);
    // we also need the identity matrix
    Operator I = Identity(FineSpace,FineSpace);
    // and we can now compute the prolongator smoother
    Operator DinvA = Dinv * FineMatrix;
    double LambdaMax = DinvA.LambdaMax();
    Operator IminusA = I - (1.333 / LambdaMax) * DinvA;
    Operator P = IminusA * Ptent;

    // define the restriction as transpose of P
    Operator R = Transpose(P);

    // use triple matrix-matrix product to define CoarseMatrix
    Operator CoarseMatrix = RAP(R,FineMatrix,P);

    // SGS will be used a smoother, and Amesos as coarse solver
    Smoother FineSolver, CoarseSolver;
    FineSolver.Reshape(FineMatrix,"SGS",MLList);
    CoarseSolver.Reshape(CoarseMatrix,"Amesos",MLList);

    // We can now construct a Preconditioner-derived object, that
    // implements the 2-level hybrid domain decomposition preconditioner.
    // Preconditioner `TwoLevelDDHybrid' can be replaced by
    // `TwoLevelDDAdditive' to define an purely additive preconditioner.
    TwoLevelDDHybrid Prec(&FineMatrix,&FineSolver,&FineSolver,
                          &CoarseSolver,&R,&P);

    // Here we define a simple Richardson method, preconditioned
    // with Prec. The exact solution is a random vector, the
    // starting solution is the zero vector. All vectors live in
    // the FineSpace.
    DoubleVector ExactSolution(FineSpace);
    DoubleVector LHS(FineSpace);
    DoubleVector RHS(FineSpace);
    DoubleVector Residual(FineSpace);
    DoubleVector PrecResidual(FineSpace);

    ExactSolution.Random();
    RHS = FineMatrix * ExactSolution;
    LHS = 0.0;
    
    double OldNorm = 1.0;

    int MaxIters = 30;
    double Tolerance = 1e-13;

    for (int i = 0 ; i < MaxIters ; ++i) {
      Residual = RHS - FineMatrix * LHS;
      Prec.Solve(Residual,PrecResidual);
      LHS = LHS + PrecResidual;
      double NewNorm = Residual.Norm2();
      if (MyPID() == 0)
        cout << "||r|| = " << Residual.Norm2() << ", "
             << "reduction = " << NewNorm / OldNorm << endl;
      if (NewNorm < Tolerance)
        break;
      OldNorm = NewNorm;
    }

  }
  catch (const char e[]) {
    cerr << "Caught exception: " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

#ifdef EPETRA_MPI
  // finalize the MLAPI workspace
  MLAPI::Finalize();
  MPI_Finalize() ;
#endif

  return 0 ;
  
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
