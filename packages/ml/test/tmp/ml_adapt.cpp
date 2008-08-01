/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_MultiLevelPreconditioner.h"

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
  // Various matrix examples are supported; please refer to the
  // Trilinos tutorial for more details.
  
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 100);
  
  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  double * x_coord = 0;
  double * y_coord = 0;
  double * z_coord = 0; // the problem is 2D, here z_coord will be NULL

  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("aggregation: damping factor", 0.0);

  // number of relaxation sweeps
  MLList.set("adaptive: max sweeps", 10);
  // number of additional null space vectors to compute
  MLList.set("adaptive: num vectors",2);

  MLList.set("adaptive: visualize", true);
  MLList.set("viz: x-coordinates", x_coord);
  MLList.set("viz: y-coordinates", y_coord);

#if 1
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, false);

  // need to allocate and fill the null space (also the
  // default one, as in this case). This vector is no longer
  // needed after a call to ComputeAdaptivePreconditioner().
  int NullSpaceSize = 1;
  vector<double> NullSpace(A->NumMyRows());
  for (int i = 0 ; i < A->NumMyRows() ; ++i)
    NullSpace[i] = 1.0;

  MLPrec->ComputeAdaptivePreconditioner(NullSpaceSize,&NullSpace[0]);
#else
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
#endif

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(1550, 1e-5);

  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(0);
  
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos");
  puts("--enable-aztecoo --enable-triutils");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */

