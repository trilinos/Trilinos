//
//  CallAmesosFactory.ccp shows how to call Amesos_Factory() with 
//  different right hand sides for each solve.  It performs three solves:
//  Solve Ax = b
//  Solve Ax1 = x
//  Solve Ax2 = x1
//
#include "Amesos_config.h"
#include "Amesos_Factory.h"
#include "Amesos_Parameter_List.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "CreateTridi.h"

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI
  // Initialize MPI

  MPI_Init(&argc,&argv);

  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm( );
#endif
  int iam = Comm.MyPID() ; 

  const int NumPoints = 10;  // Must be between 2 and 100
                              // larger numbers wil work, but the 
                              // problem is quite ill-conditioned

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.

  Epetra_Map Map(NumPoints, 0, Comm);

  //
  //  Create an empty EpetraCrsMatrix 
  Epetra_CrsMatrix A(Copy, Map, 0);

  //
  //  Populate A with a [-1,2,-1] tridiagonal matrix
  //  See CreateTridi.cpp in this directory 
  CreateTridi( A ) ; 
  

  // Create x and b vectors
  Epetra_Vector x2(Map);
  Epetra_Vector x1(Map);
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);
  Epetra_Vector residual(Map);
  Epetra_Vector temp(Map);


  Epetra_LinearProblem Problem;
  
  b.Random();

  //
  //  Solve Ax = b using KLU 
  //
  AMESOS::Parameter::List ParamList ;
  Amesos_BaseSolver* Abase ; 
  Amesos_Factory Afactory;
  Abase = Afactory.Create( AMESOS_KLU, Problem, ParamList ) ; 
  if ( Abase == 0 ) {
    cout << " AMESOS_KLU not implemented " << endl ; 
    exit(13);
  }

  Problem.SetOperator( &A );
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  //
  //  Solve Ax = b 
  //
  Problem.SetLHS( &x );
  Problem.SetRHS( &b );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //
  //  Solve Ax1 = x 
  //
  Problem.SetLHS( &x1 );
  Problem.SetRHS( &x );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //
  //  Solve Ax2 = x1
  //
  Problem.SetLHS( &x2 );
  Problem.SetRHS( &x1 );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //
  //  Compute the residual: A^3 x2 - b
  //

  A.Multiply( false, x2, temp ) ; //  temp = A x2
  A.Multiply( false, temp, x2 ) ; //  x2 = A^2 x2
  A.Multiply( false, x2, temp ) ; //  temp = A^3 x2
  residual.Update( 1.0, temp, -1.0, b, 0.0 ) ;

  double norm_residual ;
  residual.Norm2( &norm_residual ) ; 

  if (iam == 0 ) {
    cout << " norm2(A^3 x-b) = " << norm_residual << endl ; 
    //
    //  This is an ill-conditioned problem
    //
    if ( norm_residual < (1e-15)*(1.0*NumPoints*NumPoints*NumPoints*
				  NumPoints*NumPoints*NumPoints) )
      cout << " Test Passed " << endl ;
    else
      cout << " TEST FAILED " << endl ;
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ; 
#endif
  return 0;
}

