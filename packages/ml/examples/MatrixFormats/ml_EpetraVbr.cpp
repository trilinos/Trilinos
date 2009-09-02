/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

// Goal of this example if to show how to use ML with Epetra_VbrMatrix's.
//
// \author Marzio Sala, SNL 9214
// \date Last modified on 19-Jan-05

#include "ml_include.h"

// the following code cannot be compiled without these Trilinos
// packages. Note that Galeri is required in the examples only (to
// generate the linear system), not by the ML library
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"

#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

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

  // Creates the linear problem using the Galeri package.
  // Various matrix examples are supported; please refer to the
  // Galeri documentation for more details.
  // This matrix is a simple VBR matrix, constructed by replicating
  // a point-matrix on each unknown. This example is
  // useful to test the vector capabilities of ML, or to debug 
  // code for vector problems. The number of equations is here
  // hardwired as 5, but any positive number (including 1) can be
  // specified.
  //
  // NOTE: The epetra interface of ML has only limited capabilites
  // for matrices with variable block size (that is, 
  // best is if the number of equations for 
  // each block row is constant). If you are interested in variable
  // block capabilites, please contact the ML developers.
  
  int NumPDEEqns = 5;

  // build up a 9-point Laplacian in 2D. This stencil will lead to
  // "perfect" aggregates, of square shape, using almost all the ML
  // aggregation schemes.
  // The problem size (900) must be a square number. Otherwise, the user
  // can specify the number of points in the x- and y-direction, and the
  // length of the x- and y-side of the computational domain. Please
  // refer to the Trilinos tutorial for more details.
  //
  // Note also that this gallery matrix have no boundary nodes.

  int nx;
  if (argc > 1)
    nx = (int) strtol(argv[1],NULL,10);
  else
    nx = 16;
  
  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Star2D", Map, GaleriList);
  Epetra_VbrMatrix* A = CreateVbrMatrix(CrsA, NumPDEEqns);

  Epetra_Vector LHS(A->Map()); LHS.Random();
  Epetra_Vector RHS(A->Map()); RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // Construct a solver object for this problem
  AztecOO solver(Problem);

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

  // set different aggregation schemes for each level. Depending on the
  // size of your problem, the hierarchy will contain different number
  // of levels. As `Uncoupled' and `METIS' are local aggregation
  // schemes, they should be used only for the finest level. `MIS' and
  // `ParMETIS' are global aggregation schemes (meaning that the
  // aggregates can span across processes), and should be reserved for
  // coarsest levels. 
  // Note also that `Uncoupled' and `MIS' will always try to create
  // aggregates of diameter 3 (in the graph sense), while `METIS' and
  // `ParMETIS' can generate aggregates of any size.

  MLList.set("aggregation: type (level 0)", "Uncoupled");
  MLList.set("aggregation: type (level 1)", "MIS");
  
  // this is recognized by `METIS' and `ParMETIS' only
  MLList.set("aggregation: nodes per aggregate", 9);
  
  // smoother is Gauss-Seidel. Example file 
  // ml_2level_DD.cpp shows how to use
  // AZTEC's preconditioners as smoothers
  MLList.set("smoother: type","Gauss-Seidel");

  // use both pre and post smoothing. Non-symmetric problems may take
  // advantage of pre-smoothing or post-smoothing only.
  MLList.set("smoother: pre or post", "both");
  
  // solve with serial direct solver KLU
  MLList.set("coarse: type","Amesos-KLU");
  
  // create the preconditioner object and compute hierarchy
  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-5 tolerance  
  solver.Iterate(500, 1e-7);

  delete MLPrec;
  
  // compute the real residual. 

  double residual;
  LHS.Norm2(&residual);
  
  if (Comm.MyPID() == 0) 
  {
    cout << "||b-Ax||_2 = " << residual << endl;
  }

  if (residual > 1e-3)
    exit(EXIT_FAILURE);

  delete A;
  delete CrsA;
  delete Map;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
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

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-galeri");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_FAILURE);
}
#endif
