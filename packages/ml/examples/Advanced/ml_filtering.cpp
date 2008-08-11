/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
//@HEADER
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_ANASAZI)

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

#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

using namespace Teuchos;
using namespace Galeri;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  // Creates the linear problem using the Galeri package.
  // Various matrix examples are supported; please refer to the
  // Galeri documentation for more details. In this example, the grid
  // has nx x ny nodes, divided into mx x my subdomains, 
  // each assigned to a different processor.
  int nx = 8;
  int ny = 8 * Comm.NumProc();

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Recirc2D", Map, GaleriList);

  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // Construct a solver object for this problem
  AztecOO solver(Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);

  // enable filtering (also known as GGB)
  MLList.set("filtering: enable", true);
  // number of modes to be filtered
  MLList.set("filtering: eigenvalues to compute", 2);
  // number of local aggregates, the higher the better the preconditioner
  // (tough it might be expensive to apply the filtering correction)
  MLList.set("filtering: local aggregates", 1);

  // require to tune eigen-solution (here with Anasazi)
  MLList.set("eigen-analysis: tolerance", 1e-2);
  MLList.set("eigen-analysis: block-size", 1);
  MLList.set("eigen-analysis: length", 20);
  
  // create the preconditioner object and compute hierarchy
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(1550, 1e-5);

  delete MLPrec;
  
  // compute the real residual

  double residual;
  
  LHS.Norm2(&residual);

  if (Comm.MyPID() == 0) 
  {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

  if (residual > 1e-3)
    exit(EXIT_FAILURE);

  delete A;
  delete Map;

#ifdef EPETRA_MPI
  MPI_Finalize();
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
  puts("--enable-anasazi");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}

#endif 
