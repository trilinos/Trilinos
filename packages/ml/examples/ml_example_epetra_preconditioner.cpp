#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include "ml_config.h"
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

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

#include "Trilinos_Util_MatrixGallery.h"

using namespace Teuchos;

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  // Create the linear problem using the class `Trilinos_Util_MatrixGallery.'
  // Various matrix examples are supported; please refer to file
  // $TRILINOS_HOME/packages/triutils/src/Trilinos_Util_MatrixGallery.h
  // for more details.
  
  Trilinos_Util_MatrixGallery G("laplace_2d", Comm);
  G.Set("problem_size", 10000);
  
  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix * A = G.GetMatrix();
  Epetra_LinearProblem * Problem = G.GetLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  // AZTEC settings
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  MLList.setParameter("max levels",2);
  MLList.setParameter("increasing or decreasing","decreasing");

  MLList.setParameter("aggregation: type", "METIS");
  MLList.setParameter("smoother: type","aztec");
  MLList.setParameter("aggregation: nodes per aggregate", 16);
  MLList.setParameter("smoother: pre or post", "both");
  MLList.setParameter("coarse: type","Amesos_KLU");
  MLList.setParameter("coarse: max processes", 4);
  
  // Aztec smoother requires some more vectors.
  int options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE];
  AZ_defaults(options,params);
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_icc;
  MLList.setParameter("smoother: aztec options", options);
  MLList.setParameter("smoother: aztec params", params);

  // create the preconditioner object and compute hierarchy
  Epetra_ML_Preconditioner * MLPrec = new Epetra_ML_Preconditioner(*A, MLList, true);

  // verify unused parameters
  MLPrec->PrintUnused();

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  double rthresh = 1.4;
  solver.SetAztecParam(AZ_rthresh, rthresh);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  double athresh = 10.0;
  solver.SetAztecParam(AZ_athresh, athresh);
  solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

  int Niters = 500;
  solver.SetAztecOption(AZ_kspace, 160);
   
  solver.Iterate(Niters, 1e-12);

  // print out some information about the preconditioner
  //cout << MLPrec->GetOutputList();
  
  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  G.ComputeResidual(residual);
  G.ComputeDiffBetweenStartingAndExactSolutions(diff);
  
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

int main(int argc, char *argv[])
{
  puts("Please compile with --with-ml_epetra --with-ml_teuchos --with-ml_triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
