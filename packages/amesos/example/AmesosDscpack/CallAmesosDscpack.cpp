#include "Amesos_Dscpack.h"
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

  const int NumPoints = 5; // Was 1000

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.

  Epetra_Map Map(NumPoints, 0, Comm);

  cout << " CallAmesosDscpack:: NumMyPoints = " 
       << Map.NumMyPoints() 
       << " NumGlobalPoints = " 
       << Map.NumGlobalPoints() << endl ; 

  cout << " CallAmesosDscpack:: NumMyElements = " 
       << Map.NumMyElements() 
       << " NumGlobalElements = " 
       << Map.NumGlobalElements() << endl ; 


  return 0; 
  //
  //  Create an empty EpetraCrsMatrix 
  Epetra_CrsMatrix A(Copy, Map, 0);

  //
  //  Populate A with a [-1,2,-1] tridiagonal matrix
  //  See CreateTridi.cpp in this directory 
  CreateTridi( A ) ; 
  

  // Create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);
  Epetra_Vector residual(Map);
  Epetra_Vector temp(Map);


#define AXB_PROBLEM
#ifdef AXB_PROBLEM
  Epetra_LinearProblem Problem(   &A, &x, &b ) ;
#else
  Epetra_LinearProblem Problem;
#endif
  
  b.Random();

  //
  //  Solve Ax = b using DSCPACK
  //
  AMESOS::Parameter::List ParamList ;
  Amesos_Dscpack A_dscpack( Problem, ParamList ) ; 

#ifndef AXB_PROBLEM
  Problem.SetOperator( &A );
#endif
  EPETRA_CHK_ERR( A_dscpack.NumericFactorization(  ) ); 

#ifndef AXB_PROBLEM
  Problem.SetLHS( &x );
  Problem.SetRHS( &b );
#endif
  EPETRA_CHK_ERR( A_dscpack.Solve(  ) ); 

  //
  //  Compute the residual: Ax - b
  //

  A.Multiply( false, x, temp ) ; //  temp = Ax
  residual.Update( 1.0, temp, -1.0, b, 0.0 ) ;

  double norm_residual ;
  residual.Norm2( &norm_residual ) ; 

  if ( NumPoints < 10 ) {
    cout << " x = " << x << endl ; 
    cout << " b = " << b << endl ; 
  }

  if (iam == 0 ) {
    cout << " norm2(Ax-b) = " << norm_residual << endl ; 
    
    if ( norm_residual < 1e-15*NumPoints ) 
      cout << " Test Passed " << endl ;
    else
      cout << " TEST FAILED " << endl ;
  }

  return 0;
}

