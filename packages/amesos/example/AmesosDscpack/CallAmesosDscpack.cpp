//
//  CallAmesosDscpack.ccp shows how to call Amesos_Dscpack() with 
//  different right hand sides for each solve.  It performs three solves:
//  Solve Ax = b
//  Solve Ax1 = x
//  Solve Ax2 = x1
//
//
//
//  Amesos_config.h defines HAVE_AMESOS_DSCPACK if --enable-amesos-dscpack was set 
//  in the configure invocation script.
#include "Amesos_config.h"
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#include "Amesos_Parameter_List.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "CreateTridi.h"

int main(int argc, char *argv[])
{

  // Initialize MPI

  MPI_Init(&argc,&argv);

  Epetra_MpiComm Comm( MPI_COMM_WORLD );
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


#ifdef HAVE_AMESOS_DSCPACK
  Epetra_LinearProblem Problem;
  
  b.Random();

  //
  //  Solve Ax = b using DSCPACK
  //
  AMESOS::Parameter::List ParamList ;
  Amesos_Dscpack A_dscpack( Problem, ParamList ) ; 

  Problem.SetOperator( &A );
  EPETRA_CHK_ERR( A_dscpack.NumericFactorization(  ) ); 

  //
  //  Solve Ax = b 
  //
  Problem.SetLHS( &x );
  Problem.SetRHS( &b );
  EPETRA_CHK_ERR( A_dscpack.Solve(  ) ); 

  //
  //  Solve Ax1 = x 
  //
  Problem.SetLHS( &x1 );
  Problem.SetRHS( &x );
  EPETRA_CHK_ERR( A_dscpack.Solve(  ) ); 

  //
  //  Solve Ax2 = x1
  //
  Problem.SetLHS( &x2 );
  Problem.SetRHS( &x1 );
  EPETRA_CHK_ERR( A_dscpack.Solve(  ) ); 

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
#else
  cout << "This example requires DSCPACK" << endl ; 
#endif
  return 0;
}

