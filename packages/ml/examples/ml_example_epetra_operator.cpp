#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

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

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "ml_epetra_utils.h"
#include "ml_epetra_operator.h"

using namespace ML_Epetra;
using namespace Teuchos;
using namespace Trilinos_Util;

#include <iostream>

// MAIN DRIVER -- example of use of ML_Epetra::MultiLevelOperator
//
// from the command line, you may try something like that:
// $ mpirun -np 4 ./ml_example_epetra_operator.exe -problem_type=laplace_3d \
//          -problem_size=1000
//
// For more options for Trilinos_Util::CrsMatrixGallery, consult the
// Trilinos 4.0 tutorial

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  Epetra_Time Time(Comm);
  
  CommandLineParser CLP(argc,argv);
  CrsMatrixGallery G("", Comm);

  // default values for problem type and size
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "laplace_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "100" ); 

  // initialize MatrixGallery object with options specified in the shell
  Gallery.Set(CLP);
  
  // get pointer to the linear system matrix
  Epetra_CrsMatrix * A = Gallery.GetMatrix();

  // get a pointer to the map
  const Epetra_Map * Map = G.GetMap();

  // get a pointer to the linear system problem
  Epetra_LinearProblem * Problem = G.GetLinearProblem();
  
  // Construct a solver object for this problem
  AztecOO solver(*Problem);
  
  // ================= MultiLevelOperator SECTION ========================

  int nLevels = 10;
  int maxMgLevels = 6;
  ML_Set_PrintLevel(10);
  ML *ml_handle;
  
  ML_Create(&ml_handle, maxMgLevels);

  // convert to epetra matrix, put finest matrix into
  // position maxMgLevels - 1 of the hierarchy
  EpetraMatrix2MLMatrix(ml_handle, maxMgLevels-1, A);
  
  // create an Aggregate object; this will contain information
  // about the aggregation process for each level
  ML_Aggregate *agg_object;
  ML_Aggregate_Create(&agg_object);
  
  // select coarsening scheme. 
  ML_Aggregate_Set_CoarsenScheme_Uncoupled(agg_object);

  // generate the hierarchy. We decided to use decreasing ordering;
  // one can also use ML_INCREASING (in this case, you need to replace
  // maxMgLevels-1 with 0 in call to EpetraMatrix2MLMatrix())
  nLevels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, maxMgLevels-1,
						ML_DECREASING, agg_object);

  // define the ID of the coarsest level
  int coarsestLevel = maxMgLevels - nLevels;

  // set up some smoothers. Here we suppose a symmetric problem
  int nits = 1;
  double dampingFactor = 0.8;
  for(int level = maxMgLevels-1; level > coarsestLevel; level--)
    ML_Gen_Smoother_MLS(ml_handle, level, ML_BOTH, 30., 3);

  // simple coarse solver. You may want to use Amesos to access
  // to a large variety of direct solvers, serial and parallel
  ML_Gen_Smoother_SymGaussSeidel(ml_handle, coarsestLevel, ML_BOTH, 
				 nits, ML_DEFAULT);
 
  // generate the solver
  ML_Gen_Solver(ml_handle, ML_MGV, maxMgLevels-1, coarsestLevel);
 
  // create an Epetra_Operrator based on the previously created
  // hierarchy
  MultiLevelOperator MLPrec(ml_handle, Comm, *Map, *Map);

  // ========== End of MultiLevelOperator SECTION ========================
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(&MLPrec);

  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 16);

  // solve with AztecOO
  solver.Iterate(500, 1e-5);

  // compute the real residual

  double residual, diff, res2;
  G.ComputeResidual(residual);
  G.ComputeDiffBetweenStartingAndExactSolutions(diff);
  
  (G.GetExactSolution())->Norm2(&res2);
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff/res2 << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
  
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --with-ml_epetra --with-ml_teuchos --with-ml_triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */

