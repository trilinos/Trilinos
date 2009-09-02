/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
//@HEADER

// Goal of this example is to show how to use IFPACK smoothers.
//
// \author Marzio Sala, SNL 9214
//
// \date Last modified on 03-Mar-06

#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// required --enable-galeri (for the definition of the linear systems)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_IFPACK)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

using namespace Teuchos;
using namespace Galeri;

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

  // Creates the linear problem using the Galeri package. 
  int nx = 32;
  int ny = 32 * Comm.NumProc();

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Recirc2D", Map, GaleriList);

  Epetra_Vector LHS(*Map); LHS.PutScalar(0.0);
  Epetra_Vector RHS(*Map); RHS.Random();
  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // create a parameter list for ML options
  ParameterList MLList;

  // =============================== //
  // Create the ML + IFPACK smoother //
  // =============================== //
  
  ML_Epetra::SetDefaults("DD",MLList);
  MLList.set("smoother: pre or post", "post");
  MLList.set("PDE equations", 1);
  
  // fix the smoother to be IFPACK; can be set using (level X) syntax
  MLList.set("smoother: type","IFPACK");

  // now we have to specify which IFPACK preconditioner should be 
  // built. Any value that is valid for the IFPACK factory. We also need 
  // to define the overlap (>= 0).

  MLList.set("smoother: ifpack type", "ILU");
  MLList.set("smoother: ifpack overlap", 1);

  // Then, all parameters can will control the definition of the IFPACK
  // smoother are inserted in IFPACKList. In this case, we specify the fill-in
  // factor. For a list of supported parameters, please consult the IFPACK
  // documentation. For example, IFPACK preconditioner "Amesos" or
  // "Amesos stand-alone" can be used to solve with an LU 
  // factorization on each domain.
  MLList.sublist("smoother: ifpack list").set("fact: level-of-fill", 5);

  // we can now build the preconditioner...

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // ... and solve the linear system

  AztecOO Solver(Problem);
  Solver.SetPrecOperator(MLPrec);
  Solver.SetAztecOption(AZ_solver, AZ_gmres);
  Solver.SetAztecOption(AZ_kspace, 60);
  Solver.SetAztecOption(AZ_output, 5);
  Solver.Iterate(1550, 1e-5);

  // clean up and go home
  delete MLPrec;

  delete A;
  delete Map;

#ifdef HAVE_MPI
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
  puts("--enable-galeri");
  puts("--enable-aztecoo");
  puts("--enable-ifpack");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif 
