
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

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
// includes required by ML

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiLevel.h"
#include "MLAPI_EpetraPreconditioner.h"
#include "MLAPI_MATLABStream.h"
#include "MLAPI_Workspace.h"

using namespace Teuchos;
using namespace Trilinos_Util;
using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 10000);
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  AztecOO solver(*Problem);

  // =========================== begin of ML part ===========================
  
  ParameterList MLList;
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("max levels",3);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("aggregation: damping factor", 1.333);
  MLList.set("smoother: type","symmetric Gauss-Seidel");
  MLList.set("smoother: sweeps",1);
  MLList.set("smoother: damping factor",1.0);
  MLList.set("coarse: max size",32);
  MLList.set("smoother: pre or post", "both");
  MLList.set("coarse: type","Amesos-KLU");
  
  SetPrintLevel(10);
  Init();
  int size = A->NumMyRows();
  Space FineSpace(size);
  Operator AA(FineSpace,FineSpace,*A);

  MultiLevel* Cycle = new MultiLevel(AA,MLList);
  Epetra_Operator* MLAPIPrec = new EpetraPreconditioner(Comm,A->RowMatrixRowMap(),*Cycle);

  solver.SetPrecOperator(MLAPIPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 1);

  // solve with 500 iterations and 1e-12 tolerance  
  // The problem should converge as follows:
  solver.Iterate(500, 1e-5);

  delete MLAPIPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

#ifdef EPETRA_MPI
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
