#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

const bool verbose         = true;
const int OutputLevel      = 2;

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"


#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include "Teuchos_ParameterList.hpp"
#include "ml_epetra_preconditioner.h"
#include "AztecOO.h"

#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Trilinos_Util_CommandLineParser.h"


void PrintLine() 
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;
  
  return;
}

int TestMultiLevelPreconditioner(char ProblemType[],
				 Teuchos::ParameterList & MLList,
				 Epetra_LinearProblem & Problem, double & TotalErrorResidual,
				 double & TotalErrorExactSol)
{
  
  Epetra_MultiVector   * lhs     = Problem.GetLHS();
  Epetra_MultiVector   * rhs     = Problem.GetRHS();
  Epetra_RowMatrix     * A       = Problem.GetMatrix();
  
  // ======================================== //
  // create a rhs corresponding to lhs or 1's //
  // ======================================== //
  
  lhs->PutScalar(1.0);
  A->Multiply(false,*lhs,*rhs);

  lhs->PutScalar(0.0);
  
  Epetra_Time Time(A->Comm());
  
  // =================== //
  // call ML and AztecOO //
  // =================== //
  
  AztecOO solver(Problem);
  
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  
  solver.SetAztecOption(AZ_kspace, 160);
  
  solver.Iterate(1550, 1e-12);
  
  delete MLPrec;
  
  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  
  double d = 0.0, d_tot = 0.0;
  
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - 1.0) * ((*lhs)[0][i] - 1.0);
  
  A->Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if( A->Comm().MyPID() == 0 ) {
    cout << msg << "......Using " << A->Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return( solver.NumIters() );
  
}

using namespace Trilinos_Util;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumProcs = Comm.NumProc();
  
  // ===================== //
  // create linear problem //
  // ===================== //
  
  CommandLineParser CLP(argc,argv);
  CrsMatrixGallery Gallery("", Comm);
  
  // default values for problem type and size
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "recirc_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "10000" );
  
  Gallery.Set(CLP);
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  int TotalFailed = 0;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;

  int iters;
  
  // =============== //
  // default options //
  // ================ //
  if( 0 ) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    iters = TestMultiLevelPreconditioner("no defaults", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );
    cout << iters << endl;
    
    // expected iterations
    switch( NumProcs ) {
    case 1:
      if( iters != 937 ) ++TotalFailed;
      break;
    case 4:
      if( iters != 928 ) ++TotalFailed;
      break;
    case 16:
      if( iters != 1212 ) ++TotalFailed;
      break;
    default:
      if( Comm.MyPID() == 0 ) cerr << "### TEST FAILED !!!!!" << endl;
    }
  }

  // ====================== //
  // default options for SA //
  // ====================== //
  if( 1 ) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    iters = TestMultiLevelPreconditioner("SA", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    // expected iterations
    // do be fixed: iterations here may something change.
    // tolerance is to be moved to 10^10

  }

  // ====================== //
  // default options for DD //
  // ====================== //
  if( 1 ) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD",MLList);
    iters = TestMultiLevelPreconditioner("DD", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    // expected iterations
    switch( NumProcs ) {
    case 1:
      if( iters != 70 ) ++TotalFailed;
      break;
    case 4:
      if( iters != 82 ) ++TotalFailed;
      break;
    case 16:
      if( iters != 129 ) ++TotalFailed;
      break;    
    default:
      if( Comm.MyPID() == 0 ) cerr << "### TEST FAILED !!!!!" << endl;
    }
    
  }

#ifdef HAVE_ML_METIS
  // ========================================== //
  // default options for DD -- 16 aggr per proc //
  // ========================================== //
  if( 1 ) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD",MLList);
    MLList.set("aggregation: local aggregates", 16);
    iters = TestMultiLevelPreconditioner("DD", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    // expected iterations
    switch( NumProcs ) {
    case 1:
      if( iters != 65 ) ++TotalFailed;
      break;
    case 4:
      if( iters != 66 ) ++TotalFailed;
      break;
    case 16:
      if( iters != 67 ) ++TotalFailed;
      break;
    default:
      if( Comm.MyPID() == 0 ) cerr << "### TEST FAILED !!!!!" << endl;
    }
    
  }
#endif

  // ========================= //
  // default options for DD-ML //
  // ========================= //
  if( 1 ) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD-ML",MLList);
    iters = TestMultiLevelPreconditioner("DD-ML", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    // expected iterations
    switch( NumProcs ) {
    case 1:
      if( iters != 65 ) ++TotalFailed;
      break;
    case 4:
      if( iters != 77 ) ++TotalFailed;
      break;
    case 16:
      if( iters != 128 ) ++TotalFailed;
      break;
    default:
      if( Comm.MyPID() == 0 ) cerr << "### TEST FAILED !!!!!" << endl;
    }
    
  }

#if defined(HAVE_ML_METIS) && defined(HAVE_ML_PARMETIS_3x)
  // ========================= //
  // default options for DD-ML //
  // ========================= //
  {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD-ML",MLList);
    MLList.set("aggregation: nodes per aggregate (level 0)", 64);
    MLList.set("aggregation: nodes per aggregate (level 1)", 27);
    iters = TestMultiLevelPreconditioner("DD-ML", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    // expected iterations
    switch( NumProcs ) {
    case 1:
      if( iters != 51 ) ++TotalFailed;
      break;
    case 4:
      if( iters != 62 ) ++TotalFailed;
      break;
    case 16:
      if( iters != 101 ) ++TotalFailed;
      break;
    default:
      if( Comm.MyPID() == 0 ) cerr << "### TEST FAILED !!!!!" << endl;
    }
    
  }
#endif
  
  // ===================== //
  // print out total error //
  // ===================== //
  
  if( Comm.MyPID() == 0 ) {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << "......Total # of failed tests         = " << TotalFailed << endl;
    cout << endl;
 }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if( TotalFailed ) return( EXIT_FAILURE );
  else              return( EXIT_SUCCESS );

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
