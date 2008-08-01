/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
//@HEADER

// Goal of this example is to present the basic usage of
// the ML_Epetra::MultiLevelPreconditioner class.
// The example builds a simple matrix and solves the corresponding
// linear system using AztecOO and ML as a preconditioner. It finally
// checks the accuracy of the computed solution.
//
// \author Marzio Sala, ETHZ/COLAB
//
// \data Last modified on 28-Oct-05

#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example also
// requires --enable-galeri (for the definition of the linear systems)
// and --enable-aztecoo (to solve the linear system)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

// epetra objects
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
// required to build the example matrix
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
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

#ifdef ML_SCALING
   const int ntimers=4;
   enum {total, probBuild, precBuild, solve};
   ml_DblLoc timeVec[ntimers], maxTime[ntimers], minTime[ntimers];

  for (int i=0; i<ntimers; i++) timeVec[i].rank = Comm.MyPID();
  timeVec[total].value = MPI_Wtime();
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
  int nx;
  if (argc > 1)
    nx = (int) strtol(argv[1],NULL,10);
  else
    nx = 256;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

#ifdef ML_SCALING
  timeVec[probBuild].value = MPI_Wtime();
#endif
  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);
    
  // Build a linear system with trivial solution, using a random vector
  // as starting solution.
  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // As we wish to use AztecOO, we need to construct a solver object 
  // for this problem
  AztecOO solver(Problem);
#ifdef ML_SCALING
  timeVec[probBuild].value = MPI_Wtime() - timeVec[probBuild].value;
#endif

  // =========================== begin of ML part ===========================
  
#ifdef ML_SCALING
  timeVec[precBuild].value = MPI_Wtime();
#endif
  // create a parameter list for ML options
  ParameterList MLList;

  // Sets default parameters for classic smoothed aggregation. After this
  // call, MLList contains the default values for the ML parameters,
  // as required by typical smoothed aggregation for symmetric systems.
  // Other sets of parameters are available for non-symmetric systems
  // ("DD" and "DD-ML"), and for the Maxwell equations ("maxwell").
  ML_Epetra::SetDefaults("SA",MLList);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some of the parameters do not differ from their default value,
  // and they are here reported for the sake of clarity
  
  // output level, 0 being silent and 10 verbose
  MLList.set("ML output", 10);
  // maximum number of levels
  MLList.set("max levels",5);
  // set finest level to 0
  MLList.set("increasing or decreasing","increasing");

  // use Uncoupled scheme to create the aggregate
  MLList.set("aggregation: type", "Uncoupled");

  // smoother is Chebyshev. Example file 
  // `ml/examples/TwoLevelDD/ml_2level_DD.cpp' shows how to use
  // AZTEC's preconditioners as smoothers

  MLList.set("smoother: type","Chebyshev");
  MLList.set("smoother: sweeps",3);

  // use both pre and post smoothing
  MLList.set("smoother: pre or post", "both");

#ifdef HAVE_ML_AMESOS
  // solve with serial direct solver KLU
  MLList.set("coarse: type","Amesos-KLU");
#else
  // this is for testing purposes only, you should have 
  // a direct solver for the coarse problem (either Amesos, or the SuperLU/
  // SuperLU_DIST interface of ML)
  MLList.set("coarse: type","Jacobi");
#endif

  // Creates the preconditioning object. We suggest to use `new' and
  // `delete' because the destructor contains some calls to MPI (as
  // required by ML and possibly Amesos). This is an issue only if the
  // destructor is called **after** MPI_Finalize().
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // verify unused parameters on process 0 (put -1 to print on all
  // processes)
  MLPrec->PrintUnused(0);
#ifdef ML_SCALING
  timeVec[precBuild].value = MPI_Wtime() - timeVec[precBuild].value;
#endif

  // ML allows the user to cheaply recompute the preconditioner. You can
  // simply uncomment the following line:
  // 
  // MLPrec->ReComputePreconditioner();
  //
  // It is supposed that the linear system matrix has different values, but
  // **exactly** the same structure and layout. The code re-built the
  // hierarchy and re-setup the smoothers and the coarse solver using
  // already available information on the hierarchy. A particular
  // care is required to use ReComputePreconditioner() with nonzero
  // threshold.

  // =========================== end of ML part =============================
  
  // tell AztecOO to use the ML preconditioner, specify the solver 
  // and the output, then solve with 500 maximum iterations and 1e-12 
  // of tolerance (see AztecOO's user guide for more details)
  
#ifdef ML_SCALING
  timeVec[solve].value = MPI_Wtime();
#endif
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);
#ifdef ML_SCALING
  timeVec[solve].value = MPI_Wtime() - timeVec[solve].value;
#endif

  // destroy the preconditioner
  delete MLPrec;
  
  // compute the real residual

  double residual;
  LHS.Norm2(&residual);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
  }

  // for testing purposes
  if (residual > 1e-5)
    exit(EXIT_FAILURE);

  delete A;
  delete Map;

#ifdef ML_SCALING
  timeVec[total].value = MPI_Wtime() - timeVec[total].value;

  //avg
  double dupTime[ntimers],avgTime[ntimers];
  for (int i=0; i<ntimers; i++) dupTime[i] = timeVec[i].value;
  MPI_Reduce(dupTime,avgTime,ntimers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for (int i=0; i<ntimers; i++) avgTime[i] = avgTime[i]/Comm.NumProc();
  //min
  MPI_Reduce(timeVec,minTime,ntimers,MPI_DOUBLE_INT,MPI_MINLOC,0,MPI_COMM_WORLD);
  //max
  MPI_Reduce(timeVec,maxTime,ntimers,MPI_DOUBLE_INT,MPI_MAXLOC,0,MPI_COMM_WORLD);

  if (Comm.MyPID() == 0) {
    printf("timing :  max (pid)  min (pid)  avg\n");
    printf("Problem build         :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[probBuild].value,maxTime[probBuild].rank,
             minTime[probBuild].value,minTime[probBuild].rank,
             avgTime[probBuild]);
    printf("Preconditioner build  :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[precBuild].value,maxTime[precBuild].rank,
             minTime[precBuild].value,minTime[precBuild].rank,
             avgTime[precBuild]);
    printf("Solve                 :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[solve].value,maxTime[solve].rank,
             minTime[solve].value,minTime[solve].rank,
             avgTime[solve]);
    printf("Total                 :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[total].value,maxTime[total].rank,
             minTime[total].value,minTime[total].rank,
             avgTime[total]);
  }
#endif

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

#endif
