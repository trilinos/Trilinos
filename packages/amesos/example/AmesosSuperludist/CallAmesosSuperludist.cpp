// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//
//  CallAmesosSuperludist.ccp shows how to call Amesos_Factory() with 
//  Superludist, allowing the entries of A to change, but not the 
//  non-zero pattern.
//
//  Solve A x = b
//  Solve A' x1 = x
//  Solve A" x2 = x1
//  Where A' = A except that A'(0,0) = A(0,0) + 1 
//  Where A" = A' except that A"(0,0) = A'(0,0) + 1 
//
#include "Amesos_config.h"
#include "Amesos_Factory.h"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "CreateTridi.h"

#ifdef HAVE_AMESOS_SUPERLUDIST
#include "superlu_ddefs.h"
#include "supermatrix.h"
#endif


int SubTest(Epetra_Comm &Comm, Teuchos::ParameterList ParamList )
{

  int iam = Comm.MyPID() ; 
  int errors = 0 ; 
  bool verbose = false; 

  const int NumPoints = 4;  // Must be between 2 and 100 (on large matrices,
                             // the problem is quite ill-conditioned) 

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Epetra_Map Map(NumPoints, 0, Comm);

  //  Create an empty EpetraCrsMatrix 
  Epetra_CrsMatrix A(Copy, Map, 0);

  //
  //  Populate A with a [-1,2,-1] tridiagonal matrix
  //  See CreateTridi.cpp in this directory 
  CreateTridi( A ) ; 

  
  Epetra_Vector x2(Map), x1(Map), x(Map), b(Map), residual(Map), temp(Map);

  //
  //  Solve Ax = b using Amesos_KLU via the Amesos_Factory interface
  //
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos_Factory Afactory;
  //
  //  Note that Abase is created with an empty Problem, none of A, x or b
  //  have been specified at this point.  
  //  Abase = Afactory.Create( AMESOS_UMFPACK, Problem, ParamList ) ; 
  //  Abase = Afactory.Create( AMESOS_DSCPACK, Problem, ParamList ) ; 
  Abase = Afactory.Create( AMESOS_SUPERLUDIST, Problem, ParamList ) ; 
  //  Abase = Afactory.Create( AMESOS_KLU, Problem, ParamList ) ; 
  if ( Abase == 0 ) {
    cout << " AMESOS_SUPERLUDIST not implemented " << endl ; 
    exit(13);
  }

  //
  //  Factor A
  //
  Problem.SetOperator( &A );
  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  b.Random();
  b.PutScalar(1.0);
  //
  //  Solve A x = b 
  //
  Problem.SetLHS( &x );
  Problem.SetRHS( &b );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 


  //  if (verbose) cout << " x = " << x << endl ; 
  //
  int ind[1];
  double val[1];
  ind[0] = 0;
  val[0] = 1 ; 
  if ( A.MyGRID( 0 ) )
    A.SumIntoMyValues( 0, 1, val, ind ) ; 

  //  if (verbose) cout << " A' = " << A << endl ; 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  //
  //  Solve A' x1 = x 
  //
  Problem.SetLHS( &x1 );
  Problem.SetRHS( &x );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //  if (verbose) cout << " x1 = " << x1 << endl ; 

  if ( A.MyGRID( 0 ) )
    A.SumIntoMyValues( 0, 1, val, ind ) ; 

  //  if (verbose) cout << " A'' = " << A << endl ; 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  //
  //  Solve A" x2 = x1
  //
  Problem.SetLHS( &x2 );
  Problem.SetRHS( &x1 );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //  if (verbose) cout << " x2 = " << x2 << endl ; 

  //
  //  Compute the residual: A A' A" x2 - b
  //

  A.Multiply( false, x2, temp ) ; //  temp = A x2

  //  if (verbose) cout << " temp = " << temp << endl ; 

  val[0] = -val[0] ; 
  if ( A.MyGRID( 0 ) )
    A.SumIntoMyValues( 0, 1, val, ind ) ; 
  A.Multiply( false, temp, x2 ) ; //  x2 = A' A" x2



  //  if (verbose) cout << " x2 = " << x2 << endl ; 


  if ( A.MyGRID( 0 ) )
    A.SumIntoMyValues( 0, 1, val, ind ) ; 
  A.Multiply( false, x2, temp ) ; //  temp = A A' A'' x2


  //  if (verbose) cout << " temp = " << temp << endl ; 
  //  if (verbose) cout << " b = " << b << endl ; 




  residual.Update( 1.0, temp, -1.0, b, 0.0 ) ;
  //  if (verbose) cout << " residual = " << residual << endl ; 

  double norm_residual ;
  residual.Norm2( &norm_residual ) ; 

  if (iam == 0 ) {
    if (verbose) cout << " norm2(A A' A'' x-b) = " << norm_residual << endl ; 
    //
    //  This is an ill-conditioned problem
    //
    if ( norm_residual < (1e-15)*(1.0*NumPoints*NumPoints*NumPoints*
				  NumPoints*NumPoints*NumPoints) ) {
      if (verbose) cout << " Test Passed " << endl ;
    } else {
      if (verbose) cout << " TEST FAILED " << endl ;
      errors += 1 ; 
    }
  }




#define FACTOR_B
#ifdef FACTOR_B
  //
  //  Now we check to make sure that we can change the problem and 
  //  re-factorize.  
  //






  const int BNumPoints = NumPoints;  // Must be between 2 and 100 (on large matrices,
                             // the problem is quite ill-conditioned) 

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Epetra_Map BMap(BNumPoints, 0, Comm);

  //  Create an empty EpetraCrsMatrix 
  Epetra_CrsMatrix B(Copy, BMap, 0);

  //
  //  Populate A with a [-1,2,-1] tridiagonal matrix WITH -1 in the
  //  off diagonal corners.
  //  See CreateTridi.cpp in this directory 
  CreateTridiPlus( B ) ; 

  Epetra_Vector Bx2(BMap), Bx1(BMap), Bx(BMap), Bb(BMap), Bresidual(BMap), Btemp(BMap);

  //


  //
  //  Factor B
  //
  Problem.SetOperator( &B );
  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  Bb.Random();
  Bb.PutScalar(1.0);
  //
  //  Solve B x = b 
  //
  Problem.SetLHS( &Bx );
  Problem.SetRHS( &Bb );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 


  //  if (verbose) cout << " b = " << b << endl ; 
  //  if (verbose) cout << " x = " << x << endl ; 
  //
  ind[0] = 0;
  val[0] = 1 ; 
  if ( B.MyGRID( 0 ) )
    B.SumIntoMyValues( 0, 1, val, ind ) ; 

  //  if (verbose) cout << " B' = " << B << endl ; 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  //
  //  Solve B' x1 = x 
  //
  Problem.SetLHS( &Bx1 );
  Problem.SetRHS( &Bx );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //  if (verbose) cout << " x1 = " << x1 << endl ; 

  if ( B.MyGRID( 0 ) )
    B.SumIntoMyValues( 0, 1, val, ind ) ; 

  //  if (verbose) cout << " B'' = " << B << endl ; 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  //
  //  Solve B" x2 = x1
  //
  Problem.SetLHS( &Bx2 );
  Problem.SetRHS( &Bx1 );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //  if (verbose) cout << " x2 = " << x2 << endl ; 

  //
  //  Compute the residual: B B' B" x2 - b
  //

  B.Multiply( false, Bx2, Btemp ) ; //  temp = B x2

  //  if (verbose) cout << " temp = " << temp << endl ; 

  val[0] = -val[0] ; 
  if ( B.MyGRID( 0 ) )
    B.SumIntoMyValues( 0, 1, val, ind ) ; 
  B.Multiply( false, Btemp, Bx2 ) ; //  x2 = B' B" x2



  //  if (verbose) cout << " x2 = " << x2 << endl ; 


  if ( B.MyGRID( 0 ) )
    B.SumIntoMyValues( 0, 1, val, ind ) ; 
  B.Multiply( false, Bx2, Btemp ) ; //  temp = B B' B'' x2


  //  if (verbose) cout << " temp = " << temp << endl ; 
  //  if (verbose) cout << " b = " << b << endl ; 




  Bresidual.Update( 1.0, Btemp, -1.0, Bb, 0.0 ) ;
  //  if (verbose) cout << " residual = " << residual << endl ; 

  Bresidual.Norm2( &norm_residual ) ; 

  if (iam == 0 ) {
    if (verbose) cout << " norm2(B B' B'' x-b) = " << norm_residual << endl ; 
    //
    //  This is an ill-conditioned problem
    //
    if ( norm_residual < (1e-15)*(1.0*NumPoints*NumPoints*NumPoints*
				  NumPoints*NumPoints*NumPoints) ) {
      if (verbose) cout << " Test Passed " << endl ;
    } else {
      if (verbose) cout << " TEST FAILED " << endl ;
      errors += 1 ; 
    }
  }
#endif

  delete Abase;

  return errors;
}

int main(int argc, char *argv[])
{

  int errors = 0 ; 

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  int iam = Comm.MyPID() ; 
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
#if 1
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }

 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", true );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }

cout << "AAAAAA" << endl;  

 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }

 cout << "AAdddAAAA" << endl;  

 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  cout << "AAAAAAAAAaaaaa" << endl;

 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", true );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", true );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
cout << "AAAAAA" << endl;
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }

cout << "AajsdfhsAAAAA" << endl;
  
  
 {
   Teuchos::ParameterList ParamList ;
   ParamList.set( "Redistribute", false );
   ParamList.set( "AddZeroToDiag", false );
   Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.set( "ReuseSymbolic", false );
   SuperludistParams.set( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
 #endif
  

 if ( iam == 0 ) { 
   if ( errors == 0 ) 
     cout << " ALL TESTS PASSED - OK " <<  endl ; 
   else
     cout << errors << " tests failed " << endl ; 
 }


#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif
  return 0;
}
