
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
#include "ml_include.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiLevel.h"
#include "MLAPI_TwoLevelDDAdditive.h"
#include "MLAPI_TwoLevelDDHybrid.h"
#include "MLAPI_EpetraPreconditioner.h"

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

  // build the matrix using trilinos utils.
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 10000);
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  AztecOO solver(*Problem);

  // ================ B E G I N N I N G   O F   M L A P I ================
  
  // Initialize the MLAPI workspace
  Init();

  // put all the parameters required by the preconditioner in 
  // a Teuchos::ParameterList.
  ParameterList MLList;
  // maximum nunber of levels
  MLList.set("max levels",3);
  // aggregation scheme
  MLList.set("aggregation: type", "Uncoupled");
  // prolongator dampinf factor
  MLList.set("aggregation: damping factor", 0.0);
  // maximum dimension of the coarse problem
  MLList.set("coarse: max size",32);
  // smoother type (only few smoothers are available for MLAPI, 
  // including Jacobi, Gauss-Seidel, symmetric Gauss-Seidel).
  MLList.set("smoother: type","symmetric Gauss-Seidel");
  // sweeps for smoother
  MLList.set("smoother: sweeps",1);
  // damping factor for smoother
  MLList.set("smoother: damping factor",1.0);
  // `pre' will use only pre-smoother, `post' only post-smoother,
  // `both' pre- and post-smoother.
  MLList.set("smoother: pre or post", "both");
  // specify the coarse solve (only few solvers are available,
  // including Amesos-KLU)
  MLList.set("coarse: type","Amesos-KLU");
  
  ML_Set_PrintLevel(10);
  // define the space for finest-level operator
  Space FineSpace(A->NumMyRows());
  // wrap (in a light-weight mode) the Epetra_RowMatrix as MLAPI::Operator
  Operator FineMatrix(FineSpace,FineSpace,*A);

  // build the preconditioner
  MultiLevel  MLAPIPrec(FineMatrix,MLList);

  // wrap the MLAPI::Preconditioner object as an Epetra_Operator, so that
  // we can use it for AztecOO
  EpetraPreconditioner EpetraPrec(Comm,A->RowMatrixRowMap(),MLAPIPrec);

  // =========================== end of ML part =============================
  
  // inform AztecOO to use MLPrec in the preconditioning phase
  solver.SetPrecOperator(&EpetraPrec);
  
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 8);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(500, 1e-5);

  // finalize the MLAPI workspace
  Finalize();

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(0);
  
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
