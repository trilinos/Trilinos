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

#include "Amesos_ConfigDefs.h"
#include "Teuchos_ParameterList.hpp"
#include <string>
#include "Trilinos_Util_ReadTriples2Epetra.h"
#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Trilinos_Util.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"

#ifdef HAVE_AMESOS_KUNDERT
#include "KundertOO.h"
#endif
#ifdef TEST_SPOOLES
#include "SpoolesOO.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SLUD
#include "SuperludistOO.h"
#endif
#ifdef HAVE_AMESOS_SLUS
#include "Epetra_SLU.h"
#endif
#ifdef HAVE_AMESOS_SLUD2
#include "Superludist2_OO.h"
#endif
#ifdef TEST_AZTEC
#include "AztecOO.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "DscpackOO.h"
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif

#include "Amesos_TestSolver.h"
#include "CrsMatrixTranspose.h"
#include "SparseDirectTimingVars.h" 

#include <vector>
//
//  Amesos_TestSolver.cpp reads in a matrix in Harwell-Boeing format, 
//  calls one of the sparse direct solvers and computes the error 
//  and residual.  
//
//  It reads the matrix in on a single processor and can pass that
//  matrix to the solver or it can convert that matrix to a 
//  distributed matrix and pass the distributed matrix to the solver.
//
//  Amesos_TestSolver can test either A x = b or A^T x = b.
//  This can be a bit confusing because sparse direct solvers 
//  use compressed column storage - the transpose of Trilinos'
//  sparse row storage.
//
//  Matrices:
//    readA - Serial.  As read from the file.
//    transposeA - Serial.  The transpose of readA.
//    serialA - if (transpose) then transposeA else readA 
//    distributedA - readA distributed to all processes
//    passA - if ( distributed ) then distributedA else serialA
//
//
int Amesos_TestSolver( Epetra_Comm &Comm, char *matrix_file, 
		       SparseSolverType SparseSolver,
		       bool transpose, 
		       int special, AMESOS_MatrixType matrix_type ) {

  int iam = Comm.MyPID() ;




  //  int whatever;
  //  if ( iam == 0 )  cin >> whatever ; 
  //  Comm.Barrier();


  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
   
  string FileName = matrix_file ;
  int FN_Size = FileName.size() ; 
  string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );
  if ( LastFiveBytes == ".triU" ) { 
    // Call routine to read in unsymmetric Triplet matrix
    EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, false, Comm, readMap, readA, readx, 
						      readb, readxexact) );
  } else {
    if ( LastFiveBytes == ".triS" ) { 
      // Call routine to read in symmetric Triplet matrix
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, true, Comm, readMap, readA, readx, 
							readb, readxexact) );
    } else {
      if (  LastFourBytes == ".mtx" ) { 
	EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra( matrix_file, Comm, readMap, 
							       readA, readx, readb, readxexact) );
      } else {
	// Call routine to read in HB problem
	Trilinos_Util_ReadHb2Epetra( matrix_file, Comm, readMap, readA, readx, 
						     readb, readxexact) ;
      }
    }
  }

  Epetra_CrsMatrix transposeA(Copy, *readMap, 0);
  Epetra_CrsMatrix *serialA ; 

  if ( transpose ) {
    assert( CrsMatrixTranspose( readA, &transposeA ) == 0 ); 
    serialA = &transposeA ; 
  } else {
    serialA = readA ; 
  }

  Epetra_RowMatrix * passA = 0; 
  Epetra_RowMatrix * multiplyA = 0; 
  Epetra_Vector * passx = 0; 
  Epetra_Vector * passb = 0;
  Epetra_Vector * passxexact = 0;
  Epetra_Vector * passresid = 0;
  Epetra_Vector * passtmp = 0;

  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);

  Epetra_CrsMatrix A(Copy, map, 0);

  const Epetra_Map &OriginalMap = serialA->RowMatrixRowMap() ; 
  assert( OriginalMap.SameAs(*readMap) ); 
  Epetra_Export exporter(OriginalMap, map);
  Epetra_Export exporter2(OriginalMap, map);
  Epetra_Export MatrixExporter(OriginalMap, map);
  Epetra_CrsMatrix AwithDiag(Copy, map, 0);

  Epetra_Vector x(map);
  Epetra_Vector b(map);
  Epetra_Vector xexact(map);
  Epetra_Vector resid(map);
  Epetra_Vector readresid(*readMap);
  Epetra_Vector tmp(map);
  Epetra_Vector readtmp(*readMap);

  //  Epetra_Vector xcomp(map);      // X as computed by the solver
  bool distribute_matrix = ( matrix_type == AMESOS_Distributed ) ; 
  if ( distribute_matrix ) { 
    // Create Exporter to distribute read-in matrix and vectors
    //
    //  Initialize x, b and xexact to the values read in from the file
    //
    x.Export(*readx, exporter, Add);
    b.Export(*readb, exporter, Add);
    xexact.Export(*readxexact, exporter, Add);
    Comm.Barrier();
    
    A.Export(*serialA, exporter, Add);
    assert(A.FillComplete()==0);    
    
    //
    //  Add 0.0 to each diagonal entry to avoid empty diagonal entries;
    //  This is a workaround for Bug #614 
    //
    double zero = 0.0;

    AwithDiag.Export(*serialA, exporter2, Add);

    AwithDiag.SetTracebackMode(0);
    for ( int i = 0 ; i < map.NumGlobalElements(); i++ ) 
      if ( AwithDiag.LRID(i) >= 0 ) 
	AwithDiag.InsertGlobalValues( i, 1, &zero, &i ) ;
    AwithDiag.SetTracebackMode(1);

    assert(AwithDiag.FillComplete()==0);    
    
    Comm.Barrier();

    passA = &A; 
    multiplyA = &AwithDiag; 

    //    cout << " mulitplyA = " << AwithDiag << endl ; 

    passx = &x; 
    passb = &b;
    passxexact = &xexact;
    passresid = &resid;
    passtmp = &tmp;
  } else { 

    passA = serialA; 
    multiplyA = serialA; 
    passx = readx; 
    passb = readb;
    passxexact = readxexact;
    passresid = &readresid;
    passtmp = &readtmp;
  }

  Epetra_MultiVector CopyB( *passb ) ;


  double Anorm = passA->NormInf() ; 
  SparseDirectTimingVars::SS_Result.Set_Anorm(Anorm) ;

  Epetra_LinearProblem Problem(  (Epetra_RowMatrix *) passA, 
				 (Epetra_MultiVector *) passx, 
				 (Epetra_MultiVector *) passb );
  

  for ( int i = 0; i < 1+special ; i++ ) { 
    Epetra_Time TotalTime( Comm ) ; 
    
    if ( false ) { 
      //  TEST_UMFPACK is never set by configure
#ifdef TEST_UMFPACK
    } else if ( SparseSolver == UMFPACKOLD ) { 
      UmfpackOO umfpack( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
      
      umfpack.SetTrans( transpose ) ; 
      umfpack.Solve() ; 
#endif
#ifdef HAVE_AMESOS_SLUS
    } else if ( SparseSolver == SuperLU ) { 
      Epetra_SLU superluserial( &Problem ) ;
      
      superluserial.Solve() ; 
#endif
#ifdef HAVE_AMESOS_SLUD
    } else if ( SparseSolver == SuperLUdist ) {
      
	SuperludistOO superludist( Problem ) ; 
	superludist.SetTrans( transpose );  
	EPETRA_CHK_ERR( superludist.Solve( true ) ); 
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
    } else if ( SparseSolver == SUPERLUDIST ) {
	Teuchos::ParameterList ParamList ;
	Amesos_Superludist A_Superludist( Problem ) ; 

  //ParamList.set( "Redistribute", true );
  //ParamList.set( "AddZeroToDiag", true );
  Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
  ParamList.set( "MaxProcs", -3 );

	EPETRA_CHK_ERR( A_Superludist.SetParameters( ParamList ) ); 
	EPETRA_CHK_ERR( A_Superludist.SetUseTranspose( transpose ) ); 
	EPETRA_CHK_ERR( A_Superludist.SymbolicFactorization(  ) ); 
	EPETRA_CHK_ERR( A_Superludist.NumericFactorization(  ) ); 
	EPETRA_CHK_ERR( A_Superludist.Solve(  ) ); 
#endif
#ifdef HAVE_AMESOS_SLUD2
    } else if ( SparseSolver == SuperLUdist2 ) {
      
      Superludist2_OO superludist2( Problem ) ; 
      superludist2.SetTrans( transpose );  
      EPETRA_CHK_ERR( superludist2.Solve( true ) ); 
#endif
#ifdef HAVE_AMESOS_DSCPACK
    } else if ( SparseSolver == DSCPACKOLD ) {

      DscpackOO dscpack( Problem ) ; 
      EPETRA_CHK_ERR( dscpack.Solve( true ) ); 
    } else if ( SparseSolver == DSCPACK ) {
      
      Teuchos::ParameterList ParamList ;

      Amesos_Dscpack A_dscpack( Problem ) ; 
      EPETRA_CHK_ERR( A_dscpack.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_dscpack.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_dscpack.Solve(  ) ); 
#endif
#ifdef HAVE_AMESOS_MUMPS
    } else if ( SparseSolver == MUMPS ) {

	Teuchos::ParameterList ParamList ;
	Amesos_Mumps A_mumps( Problem ) ; 
	EPETRA_CHK_ERR( A_mumps.SetUseTranspose( transpose ) ); 
	EPETRA_CHK_ERR( A_mumps.SymbolicFactorization(  ) ); 
	EPETRA_CHK_ERR( A_mumps.NumericFactorization(  ) ); 
	EPETRA_CHK_ERR( A_mumps.Solve(  ) ); 

#endif
#ifdef HAVE_AMESOS_SUPERLU
    } else if ( SparseSolver == SUPERLU ) {

	Teuchos::ParameterList ParamList ;
	Amesos_Superlu A_superlu( Problem ) ; 
	EPETRA_CHK_ERR( A_superlu.SetUseTranspose( transpose ) ); 
	EPETRA_CHK_ERR( A_superlu.SymbolicFactorization(  ) ); 
	EPETRA_CHK_ERR( A_superlu.NumericFactorization(  ) ); 
	EPETRA_CHK_ERR( A_superlu.Solve(  ) ); 

#endif
#ifdef HAVE_AMESOS_SCALAPACK
    } else if ( SparseSolver == SCALAPACK ) {

	Teuchos::ParameterList ParamList ;
	Amesos_Scalapack A_scalapack( Problem ) ; 
	EPETRA_CHK_ERR( A_scalapack.SetUseTranspose( transpose ) ); 
	EPETRA_CHK_ERR( A_scalapack.SymbolicFactorization(  ) ); 
	EPETRA_CHK_ERR( A_scalapack.NumericFactorization(  ) ); 
	EPETRA_CHK_ERR( A_scalapack.Solve(  ) ); 
#endif
#ifdef HAVE_AMESOS_UMFPACK
    } else if ( SparseSolver == UMFPACK ) {

	Teuchos::ParameterList ParamList ;
	Amesos_Umfpack A_umfpack( Problem ) ; 
	EPETRA_CHK_ERR( A_umfpack.SetUseTranspose( transpose ) ); 
	EPETRA_CHK_ERR( A_umfpack.SymbolicFactorization(  ) ); 
	EPETRA_CHK_ERR( A_umfpack.NumericFactorization(  ) ); 
	EPETRA_CHK_ERR( A_umfpack.Solve(  ) ); 
#endif
#ifdef HAVE_AMESOS_KLU
    } else if ( SparseSolver == KLU ) {

	Teuchos::ParameterList ParamList ;
	Amesos_Klu A_klu( Problem ) ; 
	EPETRA_CHK_ERR( A_klu.SetUseTranspose( transpose ) ); 
	EPETRA_CHK_ERR( A_klu.SymbolicFactorization(  ) ); 
	EPETRA_CHK_ERR( A_klu.NumericFactorization(  ) ); 
	EPETRA_CHK_ERR( A_klu.Solve(  ) ); 

#endif
#ifdef TEST_SPOOLES
    } else if ( SparseSolver == SPOOLES ) { 
      SpoolesOO spooles( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
      
      spooles.SetTrans( transpose ) ; 
      spooles.Solve() ; 
#endif
#ifdef HAVE_AMESOS_KUNDERT
    } else if ( SparseSolver == KUNDERT ) { 
      KundertOO kundert( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
      
      kundert.SetTrans( transpose ) ; 
      kundert.Solve() ; 
#endif
#ifdef TEST_SPOOLESSERIAL 
    } else if ( SparseSolver == SPOOLESSERIAL ) { 
      SpoolesserialOO spoolesserial( (Epetra_RowMatrix *) passA, 
				     (Epetra_MultiVector *) passx, 
				     (Epetra_MultiVector *) passb ) ; 
      
      spoolesserial.Solve() ;
#endif
#ifndef TEST_AZTEC
    } else if ( SparseSolver == Aztec ) { 
      SparseDirectTimingVars::log_file 
	<< "Aztec Solver won't compile for me on this platform" << endl ;
      string errormsg = "Aztec Solver won't compile for me on this platform" ;
      throw( errormsg ); 
#else
    } else if ( SparseSolver == Aztec ) { 
      AztecOO aztec( (Epetra_RowMatrix *) passA, 
		     (Epetra_MultiVector *) passx, 
		     (Epetra_MultiVector *) passb ) ; 
      
      aztec.Iterate(1000, 1.0e-10 * Anorm ) ; 
#endif
    } else { 
      SparseDirectTimingVars::log_file << "Solver not implemented yet" << endl ;
      cerr << "\n\n####################  Requested solver not available on this platform ##################### ATS\n" << endl ;
      cout << " SparseSolver = " << SparseSolver << endl ; 
      cerr << " SparseSolver = " << SparseSolver << endl ; 
    }
    
    SparseDirectTimingVars::SS_Result.Set_Total_Time( TotalTime.ElapsedTime() ); 
    //    SparseDirectTimingVars::SS_Result.Set_First_Time( 0.0 ); 
    //    SparseDirectTimingVars::SS_Result.Set_Middle_Time( 0.0 ); 
    //    SparseDirectTimingVars::SS_Result.Set_Last_Time( 0.0 ); 
  }  // end for (int i=0; i<special; i++ ) 

  //
  //  Compute the error = norm(xcomp - xexact )
  //
  double error;
  passresid->Update(1.0, *passx, -1.0, *passxexact, 0.0);

  passresid->Norm2(&error);
  SparseDirectTimingVars::SS_Result.Set_Error(error) ;

  //  passxexact->Norm2(&error ) ; 
  //  passx->Norm2(&error ) ; 

  //
  //  Compute the residual = norm(Ax - b)
  //
  double residual ; 

  multiplyA->Multiply( transpose, *passx, *passtmp);
  passresid->Update(1.0, *passtmp, -1.0, *passb, 0.0); 
  //  passresid->Update(1.0, *passtmp, -1.0, CopyB, 0.0); 
  passresid->Norm2(&residual);

  SparseDirectTimingVars::SS_Result.Set_Residual(residual) ;
    
  double bnorm; 
  passb->Norm2( &bnorm ) ; 
  SparseDirectTimingVars::SS_Result.Set_Bnorm(bnorm) ;

  double xnorm; 
  passx->Norm2( &xnorm ) ; 
  SparseDirectTimingVars::SS_Result.Set_Xnorm(xnorm) ;

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  
  Comm.Barrier();

  return 0;
}
