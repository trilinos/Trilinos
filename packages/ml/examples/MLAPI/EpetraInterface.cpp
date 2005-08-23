
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
#if defined(HAVE_ML_MLAPI) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TRIUTILS)
#include "ml_config.h"
#include "ml_common.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
// required to build the example matrix
#include "Trilinos_Util_CrsMatrixGallery.h"
// required by the linear system solver
#include "AztecOO.h"

#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiLevelSA.h"
#include "MLAPI_EpetraBaseOperator.h"
#include "MLAPI_Workspace.h"

using namespace Teuchos;
using namespace MLAPI;
using namespace Trilinos_Util;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // get the epetra matrix from the Gallery
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 10000);

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  AztecOO solver(*Problem);

  Init();

  try {

    Space S(-1, A->NumMyRows(), A->RowMatrixRowMap().MyGlobalElements());
    Operator A_MLAPI(S, S, A, false);

    Teuchos::ParameterList MLList;
    MLList.set("max levels",3);
    MLList.set("increasing or decreasing","increasing");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("aggregation: damping factor", 0.0);
    MLList.set("smoother: type","symmetric Gauss-Seidel");
    MLList.set("smoother: sweeps",1);
    MLList.set("smoother: damping factor",1.0);
    MLList.set("coarse: max size",3);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type","Amesos-KLU");

    MultiLevelSA Prec_MLAPI(A_MLAPI, MLList);

    Epetra_Operator* Prec = new EpetraBaseOperator(A->RowMatrixRowMap(), Prec_MLAPI);

    solver.SetPrecOperator(Prec);
    solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
    solver.SetAztecOption(AZ_output, 32);
    solver.Iterate(500, 1e-12);

    // destroy the preconditioner
    delete Prec;

  }
  catch (const int e) {
    cerr << "Caught exception, code = " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

  Finalize();

  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);

  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  // for testing purposes
  if (residual > 1e-5)
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
  
}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("The ML API requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif // #if defined(HAVE_ML_MLAPI)
