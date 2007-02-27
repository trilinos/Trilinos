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
//#include "Trilinos_Util_ReadTriples2Epetra.h"
//#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Trilinos_Util.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef HAVE_AMESOS_SLUD
#include "SuperludistOO.h"
#endif
#ifdef HAVE_AMESOS_SLUS
#include "Epetra_SLU.h"
#endif
#ifdef HAVE_AMESOS_SLUD2
#include "Superludist2_OO.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#ifdef HAVE_AMESOS_PARDISO
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_PASTIX
#include "Amesos_Pastix.h"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
#include "Amesos_Paraklete.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef TEST_SPOOLES
#include "SpoolesOO.h"
#endif

#include "Amesos_TestSolver.h"
#include "CrsMatrixTranspose.h"
#include "SparseDirectTimingVars.h"

#include <vector>
//
//  Amesos_TestMultiSolver.cpp reads in a matrix in Harwell-Boeing format, 
//  calls one of the sparse direct solvers, using blocked right hand sides
//  and computes the error and residual.  
//
//  TestSolver ignores the Harwell-Boeing right hand sides, creating
//  random right hand sides instead.  
//
//  Amesos_TestMultiSolver can test either A x = b or A^T x = b.
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
int Amesos_TestMultiSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves, 
		      SparseSolverType SparseSolver, bool transpose,
		      int special, AMESOS_MatrixType matrix_type ) {


  int iam = Comm.MyPID() ;

  
  //  int hatever;
  //  if ( iam == 0 )  cin >> hatever ; 
  Comm.Barrier();


  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
   
  string FileName = matrix_file ;
  int FN_Size = FileName.size() ; 
  string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );
  bool NonContiguousMap = false; 

  if ( LastFiveBytes == ".triU" ) { 
    NonContiguousMap = true; 
    // Call routine to read in unsymmetric Triplet matrix
    EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, false, Comm, readMap, readA, readx, 
						      readb, readxexact, NonContiguousMap ) );
  } else {
    if ( LastFiveBytes == ".triS" ) { 
      NonContiguousMap = true; 
      // Call routine to read in symmetric Triplet matrix
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, true, Comm, 
							readMap, readA, readx, 
							readb, readxexact, NonContiguousMap ) );
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

  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);
  Epetra_Map* map_;

  if( NonContiguousMap ) {
    //
    //  map gives us NumMyElements and MyFirstElement;
    //
    int NumGlobalElements =  readMap->NumGlobalElements();
    int NumMyElements = map.NumMyElements();
    int MyFirstElement = map.MinMyGID();
    vector<int> MapMap_( NumGlobalElements );
    readMap->MyGlobalElements( &MapMap_[0] ) ;
    Comm.Broadcast( &MapMap_[0], NumGlobalElements, 0 ) ; 
    map_ = new Epetra_Map( NumGlobalElements, NumMyElements, &MapMap_[MyFirstElement], 0, Comm);
  } else {
    map_ = new Epetra_Map( map ) ; 
  }


  // Create Exporter to distribute read-in matrix and vectors
  Epetra_Export exporter(*readMap, *map_);
  Epetra_CrsMatrix A(Copy, *map_, 0);

  Epetra_RowMatrix * passA = 0; 
  Epetra_MultiVector * passx = 0; 
  Epetra_MultiVector * passb = 0;
  Epetra_MultiVector * passxexact = 0;
  Epetra_MultiVector * passresid = 0;
  Epetra_MultiVector * passtmp = 0;

  Epetra_MultiVector x(*map_,numsolves);
  Epetra_MultiVector b(*map_,numsolves);
  Epetra_MultiVector xexact(*map_,numsolves);
  Epetra_MultiVector resid(*map_,numsolves);
  Epetra_MultiVector tmp(*map_,numsolves);

  Epetra_MultiVector serialx(*readMap,numsolves);
  Epetra_MultiVector serialb(*readMap,numsolves);
  Epetra_MultiVector serialxexact(*readMap,numsolves);
  Epetra_MultiVector serialresid(*readMap,numsolves);
  Epetra_MultiVector serialtmp(*readMap,numsolves);

  bool distribute_matrix = ( matrix_type == AMESOS_Distributed ) ; 
  if ( distribute_matrix ) { 
    //
    //  Initialize x, b and xexact to the values read in from the file
    //
    
    A.Export(*serialA, exporter, Add);
    Comm.Barrier();

    assert(A.FillComplete()==0);    
    Comm.Barrier();

    passA = &A; 
    passx = &x; 
    passb = &b;
    passxexact = &xexact;
    passresid = &resid;
    passtmp = &tmp;
  } else { 
    passA = serialA; 
    passx = &serialx; 
    passb = &serialb;
    passxexact = &serialxexact;
    passresid = &serialresid;
    passtmp = &serialtmp;
  }

  passxexact->SetSeed(131) ; 
  passxexact->Random();
  passx->SetSeed(11231) ; 
  passx->Random();

  passb->PutScalar( 0.0 );
  passA->Multiply( transpose, *passxexact, *passb ) ; 

  Epetra_MultiVector CopyB( *passb ) ;

  double Anorm = passA->NormInf() ; 
  SparseDirectTimingVars::SS_Result.Set_Anorm(Anorm) ;

  Epetra_LinearProblem Problem(  (Epetra_RowMatrix *) passA, 
				 (Epetra_MultiVector *) passx, 
				 (Epetra_MultiVector *) passb );

  double max_resid = 0.0;
  for ( int i = 0 ; i < special+1 ; i++ ) { 
    
    Epetra_Time TotalTime( Comm ) ; 
    if ( false ) { 
#ifdef TEST_UMFPACK

      unused code

    } else if ( SparseSolver == UMFPACK ) { 
      UmfpackOO umfpack( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
    
      umfpack.SetTrans( transpose ) ; 
      umfpack.Solve() ; 
#endif
#ifdef TEST_SUPERLU
    } else if ( SparseSolver == SuperLU ) { 
      SuperluserialOO superluserial( (Epetra_RowMatrix *) passA, 
				     (Epetra_MultiVector *) passx, 
				     (Epetra_MultiVector *) passb ) ; 

      superluserial.SetPermc( SuperLU_permc ) ; 
      superluserial.SetTrans( transpose ) ; 
      superluserial.SetUseDGSSV( special == 0 ) ; 
      superluserial.Solve() ; 
#endif
#ifdef HAVE_AMESOS_SLUD
    } else if ( SparseSolver == SuperLUdist ) { 
      SuperludistOO superludist( Problem ) ; 
      superludist.SetTrans( transpose ) ; 
      EPETRA_CHK_ERR( superludist.Solve( true ) ) ;
#endif 
#ifdef HAVE_AMESOS_SLUD2
    } else if ( SparseSolver == SuperLUdist2 ) { 
      Superludist2_OO superludist2( Problem ) ; 
      superludist2.SetTrans( transpose ) ; 
      EPETRA_CHK_ERR( superludist2.Solve( true ) ) ;
#endif 
#ifdef TEST_SPOOLES
    } else if ( SparseSolver == SPOOLES ) { 
      SpoolesOO spooles( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
    
      spooles.SetTrans( transpose ) ; 
      spooles.Solve() ; 
#endif
#ifdef HAVE_AMESOS_DSCPACK
    } else if ( SparseSolver == DSCPACK ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Dscpack dscpack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( dscpack.SetParameters( ParamList ) ); 
    
      EPETRA_CHK_ERR( dscpack.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_UMFPACK
    } else if ( SparseSolver == UMFPACK ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Umfpack umfpack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( umfpack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( umfpack.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( umfpack.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_KLU
    } else if ( SparseSolver == KLU ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Klu klu( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( klu.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( klu.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( klu.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( klu.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( klu.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_PARAKLETE
    } else if ( SparseSolver == PARAKLETE ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Paraklete paraklete( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( paraklete.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( paraklete.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( paraklete.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( paraklete.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( paraklete.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_SLUS
    } else if ( SparseSolver == SuperLU ) { 
      Epetra_SLU superluserial( &Problem ) ; 
      EPETRA_CHK_ERR( superluserial.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( superluserial.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( superluserial.NumericFactorization(  ) ); 

      EPETRA_CHK_ERR( superluserial.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_LAPACK
    } else if ( SparseSolver == LAPACK ) { 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "MaxProcs", -3 );
      Amesos_Lapack lapack( Problem ) ; 
      EPETRA_CHK_ERR( lapack.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( lapack.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( lapack.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( lapack.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_TAUCS
    } else if ( SparseSolver == TAUCS ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Taucs taucs( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( taucs.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( taucs.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( taucs.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( taucs.NumericFactorization( ) ); 
      EPETRA_CHK_ERR( taucs.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_PARDISO
    } else if ( SparseSolver == PARDISO ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Pardiso pardiso( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( pardiso.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( pardiso.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( pardiso.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( pardiso.NumericFactorization( ) ); 
      EPETRA_CHK_ERR( pardiso.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_PASTIX
    } else if ( SparseSolver == PASTIX ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Pastix pastix( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( pastix.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( pastix.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( pastix.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( pastix.NumericFactorization( ) ); 
      EPETRA_CHK_ERR( pastix.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_PARKLETE
    } else if ( SparseSolver == PARKLETE ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Parklete parklete( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( parklete.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( parklete.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( parklete.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( parklete.NumericFactorization( ) ); 
      EPETRA_CHK_ERR( parklete.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_MUMPS
    } else if ( SparseSolver == MUMPS ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Mumps mumps( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( mumps.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( mumps.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( mumps.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( mumps.NumericFactorization( ) ); 
      EPETRA_CHK_ERR( mumps.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_SCALAPACK
    } else if ( SparseSolver == SCALAPACK ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Scalapack scalapack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( scalapack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( scalapack.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( scalapack.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( scalapack.NumericFactorization( ) ); 
      EPETRA_CHK_ERR( scalapack.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
    } else if ( SparseSolver == SUPERLUDIST ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Superludist superludist( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( superludist.SetParameters( ParamList ) ); 

      EPETRA_CHK_ERR( superludist.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( superludist.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( superludist.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( superludist.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_SUPERLU
    } else if ( SparseSolver == SUPERLU ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Superlu superlu( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( superlu.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( superlu.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( superlu.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( superlu.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( superlu.Solve( ) ); 
#endif
#ifdef TEST_SPOOLESSERIAL 
    } else if ( SparseSolver == SPOOLESSERIAL ) { 
      SpoolesserialOO spoolesserial( (Epetra_RowMatrix *) passA, 
				     (Epetra_MultiVector *) passx, 
				     (Epetra_MultiVector *) passb ) ; 
    
      spoolesserial.Solve() ;
#endif
    } else { 
      SparseDirectTimingVars::log_file << "Solver not implemented yet" << endl ;
      cerr << "\n\n####################  Requested solver not available (Or not tested with blocked RHS) on this platform #####################\n" << endl ;
    }

    SparseDirectTimingVars::SS_Result.Set_Total_Time( TotalTime.ElapsedTime() ); 
    //    SparseDirectTimingVars::SS_Result.Set_First_Time( 0.0 ); 
    //    SparseDirectTimingVars::SS_Result.Set_Middle_Time( 0.0 ); 
    //    SparseDirectTimingVars::SS_Result.Set_Last_Time( 0.0 ); 

    //
    //  Compute the error = norm(xcomp - xexact )
    //
    vector <double> error(numsolves) ; 
    double max_error = 0.0;
  
    passresid->Update(1.0, *passx, -1.0, *passxexact, 0.0);

    passresid->Norm2(&error[0]);
    for ( int i = 0 ; i< numsolves; i++ ) 
      if ( error[i] > max_error ) max_error = error[i] ; 
    SparseDirectTimingVars::SS_Result.Set_Error(max_error) ;

    //  passxexact->Norm2(&error[0] ) ; 
    //  passx->Norm2(&error ) ; 

    //
    //  Compute the residual = norm(Ax - b)
    //
    vector <double> residual(numsolves) ; 
  
    passtmp->PutScalar(0.0);
    passA->Multiply( transpose, *passx, *passtmp);
    passresid->Update(1.0, *passtmp, -1.0, *passb, 0.0); 
    //    passresid->Update(1.0, *passtmp, -1.0, CopyB, 0.0); 
    passresid->Norm2(&residual[0]);

    for ( int i = 0 ; i< numsolves; i++ ) 
      if ( residual[i] > max_resid ) max_resid = residual[i] ; 


    SparseDirectTimingVars::SS_Result.Set_Residual(max_resid) ;
    
    vector <double> bnorm(numsolves); 
    passb->Norm2( &bnorm[0] ) ; 
    SparseDirectTimingVars::SS_Result.Set_Bnorm(bnorm[0]) ;

    vector <double> xnorm(numsolves); 
    passx->Norm2( &xnorm[0] ) ; 
    SparseDirectTimingVars::SS_Result.Set_Xnorm(xnorm[0]) ;


    if ( false && iam == 0 ) { 

      cout << " Amesos_TestMutliSolver.cpp " << endl ; 
      for ( int i = 0 ; i< numsolves && i < 10 ; i++ ) {
	cout << "i=" << i 
	     << " error = " << error[i] 
	     << " xnorm = " << xnorm[i] 
	     << " residual = " << residual[i] 
	     << " bnorm = " << bnorm[i] 
	     << endl ; 
      
      }
    
      cout << endl << " max_resid = " << max_resid ; 
      cout << " max_error = " << max_error << endl ; 
      cout << " Get_residual() again = " << SparseDirectTimingVars::SS_Result.Get_Residual() << endl ;

    }
  }
  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  delete map_;
  
  Comm.Barrier();

return 0 ;
}
