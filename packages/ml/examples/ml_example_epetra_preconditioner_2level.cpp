#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include "ml_config.h"
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_epetra_preconditioner.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Teuchos;
using namespace Trilinos_Util;

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  // Create the linear problem using the class `Trilinos_Util::CrsMatrixGallery.'
  // Various matrix examples are supported; please refer to file
  // $TRILINOS_HOME/packages/triutils/src/Trilinos_Util_CrsMatrixGallery.h
  // for more details.
  
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 10000);
  
  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("DD",MLLIst);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some parameters are reported here to better explain the process
  // even if they are as defaults
  MLList.set("aggregation: type", "METIS");
  MLList.set("smoother: type","Aztec");
  // put 16 nodes on each aggregate. This number can be too small
  // for large problems. In this case, either augment ir, or increase
  // the number of levels
  MLList.set("aggregation: nodes per aggregate", 16);
  // problem is symmetric, pre and post smoothing
  MLList.set("smoother: pre or post", "both");
  // use a serial solver, included in Amesos
  MLList.set("coarse: type","Amesos_KLU");
  
  // now we need to define the solvers on each subdomain (== processor).
  // Here we use an incomplete Cholesky factorization, with no fill-in
  // and no overlap. To that aim, we use Aztec's preconditioning function.
  // Aztec requires two more vectors. Note: the following options and params
  // will be used ONLY for the smoother, and will NOT affect the Aztec solver
  int options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE];
  AZ_defaults(options,params);
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_icc;
  // stick pointers in the parameters' list. NOTE THAT THE VECTORS ARE NOT
  // COPIED! Only the pointer is copied, so do not destroy options and params
  // before the end of the linear system solution!
  MLList.set("smoother: Aztec options", options);
  MLList.set("smoother: Aztec params", params);

  // create the preconditioner object and compute hierarchy
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // verify unused parameters, print unused parameters on process 0
  MLPrec->PrintUnused(0);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================

  // use CG, and estimate the condition number.
  // A more classical option like AZ_cg is fine too.  
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 as tolerance on the
  // relative residual  
  solver.Iterate(500, 1e-12);

  // delete the preconditioner. Do it BEFORE MPI_Finalize
  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
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
  puts("Please configure ML with --with-ml_epetra --with-ml_teuchos --with-ml_triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
