/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
//@HEADER

// Goal of this example is to show the (limited) ML analysis capabilities.
//
// \author Marzio Sala, SNL 9214
//
// \date Last modified on 14-Jun-05

#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// required --enable-galeri (for the definition of the linear systems)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) 

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
  // Several matrix examples are supported; please refer to the
  // Galeri documentation for more details.
  // Most of the examples using the ML_Epetra::MultiLevelPreconditioner
  // class are based on Epetra_CrsMatrix. Example
  // `ml_EpetraVbr.cpp' shows how to define a Epetra_VbrMatrix.
  
  // `Laplace2D' is a symmetric matrix; an example of non-symmetric
  // matrices is `Recirc2D' (advection-diffusion in a box, with
  // recirculating flow). The grid has nx x ny nodes, divided into
  // mx x my subdomains, each assigned to a different processor.
  int nx = 28;
  int ny = 28 * Comm.NumProc();

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);

  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  
  // use Uncoupled scheme to create the aggregate
  MLList.set("aggregation: type", "Uncoupled");
  
  // fix the smoother
  MLList.set("smoother: type","symmetric Gauss-Seidel");

  // create the preconditioning object.
  // Note that users need to set "viz: enable" == true in order to
  // visualize!

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // for 2D Cartesian grid, you can print the stencil of your operator
  // using this simple function.
  
  MLPrec->PrintStencil2D(nx, ny);

  // ================================================= //
  // A N A L Y S I S   O F   T H E   H I E R A R C H Y //
  // ================================================= //

  // Method AnalyzeHierarchy() can be used to validate an
  // already built hierarchy.
  // - `true' means perform a "cheap" analysis of each level's matrix
  // - Then, each level's smoothers and the complete cycle are 
  //   applied to solve the problem
  //     A e = 0
  //   with a random initial solution, to get a sense of the effectiveness
  //   of the smoothers and the cycle itself. The parameters are:
  //   * NumPreCycles and NumPostCycles specify the number of post
  //     and pre smoother applications;
  //   * NumMLCycles specifies the number of applications of the 
  //     complete cycle.

  int NumPreCycles = 5;
  int NumPostCycles = 1;
  int NumMLCycles = 10;
  MLPrec->AnalyzeHierarchy(true, NumPreCycles, NumPostCycles, NumMLCycles);

  // ================================================= //
  // A N A L Y S I S   O F   T H E   S M O O T H E R S //
  // ================================================= //

  // Method TestSmoothers() can be used to analyze different smoothers
  // on a given problem. The cycle is built following the parameters
  // specified in MLTestList.
  // Please refer to the user's guide for more details about the following
  // parameters. Not all smoothers are supported by this testing.

  Teuchos::ParameterList MLTestList;
  ML_Epetra::SetDefaults("SA",MLTestList);
  
  MLPrec->TestSmoothers(MLTestList);

  // =============== //
  // end of analysis //
  // =============== //

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

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif 
