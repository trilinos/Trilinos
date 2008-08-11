
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

// Goal of this example is to present how to define domain decomposition
// preconditioner using class ML_Epetra::MultiLevelPreconditioner.
//
// \author Marzio Sala, SNL 9214
// \date Last modified on 18-Jan-05

#include "ml_include.h"

// configure options:
// ------------------
// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), required Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// requires --enable-galeri (for the definition of the linear systems)
// and --enable-aztecoo (for the solution of the linear system).

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

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
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
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

  Epetra_Time Time(Comm);

  // Create the linear problem using the Galeri package. 
  // Here, we are using a symmetric matrix; non-symmetric matrices
  // can also be defined. Please refer to the Galeri documentation 
  // for more details.
  // Here the grid has size nx x ny x nz, and it is subdivided into
  // mx x my x mz subdomains (== processors).

  ParameterList GaleriList;
  GaleriList.set("nx", 8);
  GaleriList.set("ny", 8);
  GaleriList.set("nz", 8 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", 1);
  GaleriList.set("mz", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian3D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace3D", Map, GaleriList);

  // Construct the linear system with trivial solution
  
  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // Construct a solver object for this problem
  
  AztecOO solver(Problem);

  // =========================== begin of ML part ===========================
  
  // create an empty parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation with heavy smoothers
  // (of domain decomposition type, i.e. one-level Schwarz with incomplete
  // factorizations on each subdomain/process)
  // We need to define the solvers on each subdomain (== processor).
  // Here we use an incomplete LU factorization, with no fill-in
  // and no overlap. To that aim, we use Aztec's preconditioning function.
  // Aztec requires two more vectors. Note: the following options and params
  // will be used ONLY for the smoother, and will NOT affect the Aztec solver
  // NOTE: to use exact solvers change to AZ_lu (requires AztecOO configured
  // with option--enable-aztecoo-azlu), of use IFPACK smoothers 
  // (requires Trilinos to be built with options--enable-ifpack --enable-amesos)
  
  int options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE];
  AZ_defaults(options,params);
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_ilu;
  options[AZ_graph_fill] = 0;
  options[AZ_overlap] = 0;
  
  // SetDefaults() will call AZ_defaults(options,params), and will also set the
  // preconditioner as `AZ_dom_decomp'. 
  // NOTE THAT THE `options' AND `params' VECTORS ARE NOT COPIED into
  // the list, only the pointers is stored, so do not delete options 
  // and params before the end of the linear system solution!
  // Alternatively, you can also call SetDefaults() without passing 
  // `options' and `params.' This way, the code will allocate a int 
  // and a double vector, that must be freed by the user.
  // `DD' means to set default values for domain decomposition
  // preconditioners
  
  ML_Epetra::SetDefaults("DD",MLList,options,params);
  
  // Overwrite some parameters. Please refer to the user's guide
  // for more information
  // Some parameters are reported here to better explain the process
  // even if they are as defaults. 
  // NOTE: To use `METIS' as aggregation scheme, you need to configure
  // ML with the option --with-ml_metis. Otherwise, the code will
  // creates aggregates containing all the local nodes (that is,
  // the dimension of the coarse problem will be equal to the
  // number of processors)
 
  MLList.set("aggregation: type", "METIS");
  MLList.set("smoother: type","Aztec");
  
  // Put 64 nodes on each aggregate. This number can be too small
  // for large problems. In this case, either augment this value, or increase
  // the number of levels. Also, use only presmoothing, and KLU as
  // coarse solver (KLU is enabled by default with Amesos)
  
  MLList.set("aggregation: nodes per aggregate", 64);
  MLList.set("smoother: pre or post", "pre");
  MLList.set("coarse: type","Amesos-KLU");
  
  // Create the preconditioning object. We suggest to use `new' and
  // `delete' because the destructor contains some calls to MPI (as
  // required by ML and possibly Amesos). This is an issue only if the
  // destructor is called **after** MPI_Finalize().
 
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================

  // Instruct AztecOO to use GMRES with estimation of the condition
  // number. Also, requires output every 32 iterations
  // Then, solve with 500 iterations and 1e-12 as tolerance on the
  // relative residual  
  
  solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  // delete the preconditioner. Do it BEFORE MPI_Finalize
  delete MLPrec;
  
  // compute the real residual

  double residual;
  LHS.Norm2(&residual);
  
  if (Comm.MyPID() == 0)
  {
    cout << "||x_exact - x||_2 = " << residual << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

  delete A;
  delete Map;

  if (residual > 1e-5)
    exit(EXIT_FAILURE);

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
  
  exit(EXIT_SUCCESS);
}
#endif 
