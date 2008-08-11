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
// required to build the example matrix
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
#include "ml_MultiLevelPreconditioner.h"
#include "Galeri_ReadHB.h"

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

  // simple reading of parameters:
  // - first is nx
  // - second is problem type
  // - third is preconditioner
  
  int nx, ny;
  nx = (int) strtol(argv[1],NULL,10);
  ny = nx;

  string ProblemType = argv[2];
  cout << "PROBLEM TYPE = " << ProblemType << endl;

  string what = argv[3];
  cout << "CHOICE = " << what << endl;

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
  GaleriList.set("alpha", 3.14152 / 2);

  Epetra_Map* Map;
  Epetra_CrsMatrix* A;
  if (false)
  {
    Map = CreateMap("Cartesian2D", Comm, GaleriList);
    A = CreateCrsMatrix(ProblemType, Map, GaleriList);
  }
  else
  {
    Epetra_Vector* readx = 0,* readb = 0,* readxexact = 0;
    Galeri::ReadHB(ProblemType.c_str(), Comm, Map, A, readx, readb, readxexact);
  }
    
  // Build a linear system with trivial solution, using a random vector
  // as starting solution.
  Epetra_Vector LHS_exact(*Map); LHS_exact.Random();
  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); 
  A->Multiply(false, LHS_exact, RHS);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // As we wish to use AztecOO, we need to construct a solver object 
  // for this problem
  AztecOO solver(Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);
  
  // output level, 0 being silent and 10 verbose
  MLList.set("ML output", 10);
  // maximum number of levels
  MLList.set("max levels",5);
  // set finest level to 0
  MLList.set("increasing or decreasing","increasing");

  // use Uncoupled scheme to create the aggregate
  MLList.set("aggregation: type", "Uncoupled");

  // smoother is symmetric Gauss-Seidel. Example file 
  // `ml/examples/TwoLevelDD/ml_2level_DD.cpp' shows how to use
  // AZTEC's preconditioners as smoothers
  MLList.set("smoother: type", "symmetric Gauss-Seidel");
  MLList.set("smoother: sweeps", 2);
  MLList.set("smoother: damping factor", 0.67);

  // use both pre and post smoothing
  MLList.set("smoother: pre or post", "both");
  MLList.set("PDE equations", 1);

  if (what == "SA")
  {
    MLList.set("aggregation: damping factor", 1.333);
  }
  else if (what == "NSA")
  {
    MLList.set("aggregation: damping factor", 0.0);
  }
  else if (what == "NSR")
  {
    MLList.set("aggregation: use tentative restriction", true);
  }
  else if (what == "1")
  {
    MLList.set("energy minimization: enable", true);
    MLList.set("energy minimization: type", 1);
  }
  else if (what == "2")
  {
    MLList.set("energy minimization: enable", true);
    MLList.set("energy minimization: type", 2);
  }
  else if (what == "3")
  {
    MLList.set("energy minimization: enable", true);
    MLList.set("energy minimization: type", 3);
  }
  else if (what == "AAT")
  {
    MLList.set("aggregation: symmetrize matrix", true);
  }
  else
    exit(-1);
               
  MLList.set("coarse: type","Amesos-KLU");

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // =========================== end of ML part =============================
  
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_kspace, 60);
  solver.SetAztecOption(AZ_output, 5);
  solver.Iterate(1550, 1e-5);

  LHS_exact.Update(1.0, LHS, -1.0);
  double norm;
  LHS_exact.Norm2(&norm);

  cout << "||X - X_exact||_2 = " << norm << endl;

  // destroy the preconditioner
  delete MLPrec;
  
  delete A;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
