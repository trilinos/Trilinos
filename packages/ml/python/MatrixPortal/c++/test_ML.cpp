#include "ml_include.h"
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
#include "AztecOO.h"
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
  int nx;
  if (argc > 1)
    nx = (int) strtol(argv[1],NULL,10);
  else
    nx = 8;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

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

  // =========================== begin of ML part ===========================
  
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
  MLList.set("output", 10);
  // maximum number of levels
  MLList.set("max levels",5);
  // set finest level to 0
  MLList.set("increasing or decreasing","increasing");

  // use Uncoupled scheme to create the aggregate
  MLList.set("aggregation: type", "Uncoupled");

  // smoother is symmetric Gauss-Seidel. Example file 
  // `ml/examples/TwoLevelDD/ml_2level_DD.cpp' shows how to use
  // AZTEC's preconditioners as smoothers
  MLList.set("smoother: type","symmetric Gauss-Seidel");

  // use both pre and post smoothing
  MLList.set("smoother: pre or post", "both");

  // solve with serial direct solver KLU
  MLList.set("coarse: type","Amesos-KLU");

  // Creates the preconditioning object. We suggest to use `new' and
  // `delete' because the destructor contains some calls to MPI (as
  // required by ML and possibly Amesos). This is an issue only if the
  // destructor is called **after** MPI_Finalize().
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // =========================== end of ML part =============================
  
  // tell AztecOO to use the ML preconditioner, specify the solver 
  // and the output, then solve with 500 maximum iterations and 1e-12 
  // of tolerance (see AztecOO's user guide for more details)
  
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

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

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
