
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

// Goal of this example is to present the visualization capabilities of
// ML. Using ML, the user can visualize the aggregates for all levels.
// This requires, as additional input, the coordinates of the fine-grid
// nodes. The output file is simple collection of 2D or 3D points,
// each of them containing the (double) value of the aggregate it belongs to.
// A freely-downloadable software, called XD3D, can for example
// be used to visualize the aggregates. ML can also visualize the effect
// of smoothers and the entire ML cycle on random vectors; see the
// `visualization' section of this example.
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
#include "Galeri_Utils.h"

#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Create the linear problem using the Galeri package.
  
  int NumPDEEqns = 5;

  Teuchos::ParameterList GaleriList;
  int nx = 32;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Laplace2D", Map, GaleriList);
  Epetra_VbrMatrix* A = CreateVbrMatrix(CrsA, NumPDEEqns);

  Epetra_Vector LHS(A->Map()); LHS.Random();
  Epetra_Vector RHS(A->Map()); RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  AztecOO solver(Problem);

  // =========================== definition of coordinates =================
  
  // use the following Galeri function to get the
  // coordinates for a Cartesian grid. Note however that the
  // visualization capabilites of Trilinos accept non-structured grid as
  // well. Visualization and statistics occurs just after the ML
  // preconditioner has been build.

  Epetra_MultiVector* Coord = CreateCartesianCoordinates("2D", &(A->Map()),
                                                         GaleriList);
  double* x_coord = (*Coord)[0];
  double* y_coord = (*Coord)[1];
  
  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;
  int *options    = new int[AZ_OPTIONS_SIZE];
  double *params  = new double[AZ_PARAMS_SIZE];

  // set defaults
  ML_Epetra::SetDefaults("SA",MLList, options, params);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some of the parameters do not differ from their default value,
  // and they are here reported for the sake of clarity
  
  // maximum number of levels
  MLList.set("max levels",3);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("smoother: type", "symmetric Gauss-Seidel");

  // aggregation scheme set to Uncoupled. Note that the aggregates
  // created by MIS can be visualized for serial runs only, while 
  // Uncoupled, METIS for both serial and parallel runs.
  MLList.set("aggregation: type", "Uncoupled");

  // ======================== //
  // visualization parameters //
  // ======================== //
  // 
  // - set "viz: enable" to `false' to disable visualization and
  //   statistics.
  // - set "x-coordinates" to the pointer of x-coor
  // - set "viz: equation to plot" to the number of equation to 
  //   be plotted (for vector problems only). Default is -1 (that is,
  //   plot all the equations)
  // - set "viz: print starting solution" to print on file 
  //   the starting solution vector, that was used for pre-
  //   and post-smoothing, and for the cycle. This may help to
  //   understand whether the smoothed solution is "smooth" 
  //   or not.
  //
  // NOTE: visualization occurs *after* the creation of the ML preconditioner,
  // by calling VisualizeAggregates(), VisualizeSmoothers(), and
  // VisualizeCycle(). However, the user *must* enable visualization 
  // *before* creating the ML object. This is because ML must store some
  // additional information about the aggregates.
  // 
  // NOTE: the options above work only for "viz: output format" == "xyz"
  // (default value) or "viz: output format" == "vtk".
  // If "viz: output format" == "dx", the user
  // can only plot the aggregates. 

  MLList.set("viz: output format", "vtk");
  MLList.set("viz: enable", true);
  MLList.set("x-coordinates", x_coord);
  MLList.set("y-coordinates", y_coord);
  MLList.set("z-coordinates", (double *)0);
  MLList.set("viz: print starting solution", true);

  // =============================== //
  // end of visualization parameters //
  // =============================== //

  // create the preconditioner object and compute hierarchy

  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // ============= //
  // visualization //
  // ============= //

  // 1.- print out the shape of the aggregates, plus some
  //     statistics
  // 2.- print out the effect of presmoother and postsmoother
  //     on a random vector. Input integer number represent 
  //     the number of applications of presmoother and postmsoother,
  //     respectively
  // 3.- print out the effect of the ML cycle on a random vector.
  //     The integer parameter represents the number of cycles.
  // Below, `5' and `1' refers to the number of pre-smoother and
  // post-smoother applications. `10' refers to the number of ML
  // cycle applications. In both cases, smoothers and ML cycle are
  // applied to a random vector.

  MLPrec->VisualizeAggregates();
  MLPrec->VisualizeSmoothers(5,1);
  MLPrec->VisualizeCycle(10);

  // ==================== //
  // end of visualization //
  // ==================== //

  // destroy the preconditioner
  delete MLPrec;
  
  delete [] options;
  delete [] params;
  
  delete A;
  delete Coord;
  delete Map;

#ifdef HAVE_MPI
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

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) */
