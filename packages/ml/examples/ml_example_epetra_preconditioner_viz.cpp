#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

#include "ml_include.h"

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
  // Various matrix examples are supported; please refer to the
  // Trilinos tutorial for more details.
  // The matrix here is a VBR matrix, but the code works with any matrix
  // as well. Just note that the coordinates refer to the block rows: if
  // the matrix size if 10000*NumPDEEqns, the coordinate vectors are
  // allocated of size 10000. 
  
  int NumPDEEqns = 5;

  VbrMatrixGallery Gallery("laplace_2d_9pt", Comm);
  Gallery.Set("problem_size", 10000);

  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetVbrMatrix(NumPDEEqns);
  Epetra_LinearProblem * Problem = Gallery.GetVbrLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);
  // get coordinates of points (for Cartesian grid)

  // =========================== definition of coordinates =================
  
  double * x_coord;
  double * y_coord;
  double * z_coord; // the problem is 2D, here z_coord will be NULL
  
  // use the following triutils matrix gallery function to get the
  // coordinates for a Cartesian grid. Note however that the
  // visualization capabilites of Trilinos accept non-structured grid as
  // well. Visualization and statistics occurs just after the ML
  // preconditioner has been build.

  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some of the parameters do not differ from their default value,
  // and they are here reported for the sake of clarity
  
  // maximum number of levels
  MLList.set("max levels",5);
  MLList.set("increasing or decreasing","increasing");

  // aggregation scheme set to Uncoupled. Note that MIS can be
  // visualized only for serial runs, while Uncoupled, METIS and
  // ParMETIS for serial and parallel runs.
  MLList.set("aggregation: type", "Uncoupled");
  
  // ============== visualization with OpenDX. ==================
  // - set "viz: enable" to false to disable visualization and
  //   statistics.
  // - "viz: dimensions" is the number of dimensions, from 1 to 3
  // - "viz: z-coordinates" can be used for 3D problems.
  MLList.set("viz: enable", true);
  MLList.set("viz: dimensions", 2);
  MLList.set("viz: x-coordinates", x_coord);
  MLList.set("viz: y-coordinates", y_coord);
  // ============== end of visualization parameters =============

  // create the preconditioner object and compute hierarchy
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(500, 1e-12);

  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidualVbr(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutionsVbr(diff);
  
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
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
