#include "Amesos_ConfigDefs.h"

#include "Trilinos_Util.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef TEST_KUNDERT
#include "KundertOO.h"
#endif
#ifdef HAVE_AMESOS_SLUD
#include "SuperludistOO.h"
#endif
#ifdef HAVE_AMESOS_SLUD2
#include "Superludist2_OO.h"
#endif
#ifdef TEST_SPOOLES
#include "SpoolesOO.h"
#endif
#ifdef TEST_AZTEC
#include "AztecOO.h"
#endif
#if 0
#include "UmfpackOO.h"
#include "SpoolesserialOO.h"
#include "SuperluserialOO.h"
#endif
#include "SparseSolverResult.h"
#include "Amesos_TestSolver.h"
#include "CrsMatrixTranspose.h"
#include "SparseDirectTimingVars.h"

#include <vector>
//
//  TestMrhsSolver.cpp reads in a matrix in Harwell-Boeing format, 
//  calls one of the sparse direct solvers, using multiple right hand sides
//  (one per solve) and computes the error and residual.  
//
//  TestSolver ignores the Harwell-Boeing right hand sides, creating
//  random right hand sides instead.  
//
//  TestMrhsSolver can test either A x = b or A^T x = b.
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
int Amesos_TestMrhsSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves, 
		     SparseSolverType SparseSolver, bool transpose, 
		     int special ) {


  int iam = Comm.MyPID() ;

  //  int hatever;
  //  if ( iam == 0 )  cin >> hatever ; 
  Comm.Barrier();


  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
   
  // Call routine to read in HB problem
  Trilinos_Util_ReadHb2Epetra( matrix_file, Comm, readMap, readA, readx, 
			       readb, readxexact);


  Epetra_CrsMatrix transposeA(Copy, *readMap, 0);
  Epetra_CrsMatrix *serialA ; 

  if ( transpose ) {
    assert( CrsMatrixTranspose( readA, &transposeA ) == 0 ); 
    serialA = &transposeA ; 
  } else {
    serialA = readA ; 
  }

  


  Epetra_RowMatrix * passA = 0; 
  Epetra_MultiVector * passx = 0; 
  Epetra_MultiVector * passb = 0;
  Epetra_MultiVector * passxexact = 0;
  Epetra_MultiVector * passtemp = 0;
  Epetra_MultiVector * passtmp = 0;

  Epetra_MultiVector x(*readMap,numsolves);
  Epetra_MultiVector b(*readMap,numsolves);
  Epetra_MultiVector xexact(*readMap,numsolves);
  Epetra_MultiVector resid(*readMap,numsolves);
  Epetra_MultiVector readresid(*readMap,numsolves);
  Epetra_MultiVector tmp(*readMap,numsolves);
  Epetra_MultiVector readtmp(*readMap,numsolves);

  passA = serialA; 
  passx = &x; 
  passb = &b;
  passxexact = &xexact;
  passtemp = &readresid;
  passtmp = &readtmp;

  passxexact->SetSeed(1.31) ; 
  passxexact->Random();
  passx->SetSeed(11.231) ; 
  passx->Random();

  serialA->Multiply( false, *passxexact, *passb ) ; 

  double Anorm = passA->NormInf() ; 
  SparseDirectTimingVars::SS_Result.Set_Anorm(Anorm) ;

  Epetra_LinearProblem Problem(  (Epetra_RowMatrix *) passA, 
				 (Epetra_MultiVector *) passx, 
				 (Epetra_MultiVector *) passb );

  Epetra_Time TotalTime( Comm ) ; 
  //  switch( SparseSolver ) { 
  //  case UMFPACK:
  if ( false ) { 
#ifdef TEST_UMFPACK
  } else if ( SparseSolver == UMFPACK ) { 
    UmfpackOO umfpack( (Epetra_RowMatrix *) passA, 
		       (Epetra_MultiVector *) passx, 
		       (Epetra_MultiVector *) passb ) ; 
    
    umfpack.SetTrans( transpose ) ; 
    umfpack.Solve() ; 
#endif
#ifdef TEST_SUPERLU
  } else if ( SparseSolver == SuperLU ) { 
    SuperluserialOO superluserial ; 
    superluserial.SetUserMatrix( (Epetra_RowMatrix *) passA) ; 

    superluserial.SetPermc( SuperLU_permc ) ; 
    superluserial.SetTrans( transpose ) ; 
    superluserial.SetUseDGSSV( special == 0 ) ; 

    for ( int i= 0 ; i < numsolves ; i++ ) { 
      //    set up to sovle A X[:,i] = B[:,i]
      Epetra_Vector *passb_i = (*passb)(i) ;
      Epetra_Vector *passx_i = (*passx)(i) ;
      superluserial.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
      superluserial.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
      //      superluserial.SetRHS( (Epetra_MultiVector *) passb_i ; 
      superluserial.Solve() ; 
      if ( i == 0 ) 
	SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
      if ( i < numsolves-1 ) 
	SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
      else
	SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 

    }
#endif
#ifdef HAVE_AMESOS_SLUD
  } else if ( SparseSolver == SuperLUdist ) { 
    SuperludistOO superludist( Problem ) ; 
    superludist.SetTrans( transpose ) ; 

    bool factor = true; 
    for ( int i= 0 ; i < numsolves ; i++ ) { 
      //    set up to sovle A X[:,i] = B[:,i]
      Epetra_Vector *passb_i = (*passb)(i) ;
      Epetra_Vector *passx_i = (*passx)(i) ;
      Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
      Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
      EPETRA_CHK_ERR( superludist.Solve( factor ) ); 
      factor = false; 
      if ( i == 0 ) 
	SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
      if ( i < numsolves-1 ) 
	SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
      else
	SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 

    }
#endif
#ifdef HAVE_AMESOS_SLUD2
  } else if ( SparseSolver == SuperLUdist2 ) { 
    Superludist2_OO superludist2( Problem ) ; 
    superludist2.SetTrans( transpose ) ; 

    bool factor = true; 
    for ( int i= 0 ; i < numsolves ; i++ ) { 
      //    set up to sovle A X[:,i] = B[:,i]
      Epetra_Vector *passb_i = (*passb)(i) ;
      Epetra_Vector *passx_i = (*passx)(i) ;
      Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
      Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
      EPETRA_CHK_ERR( superludist2.Solve( factor ) ); 
      factor = false; 
      if ( i == 0 ) 
	SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
      if ( i < numsolves-1 ) 
	SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
      else
	SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 

    }
#endif
#ifdef TEST_SPOOLES
  } else if ( SparseSolver == SPOOLES ) { 
    SpoolesOO spooles( (Epetra_RowMatrix *) passA, 
		       (Epetra_MultiVector *) passx, 
		       (Epetra_MultiVector *) passb ) ; 
    
    spooles.SetTrans( transpose ) ; 
    spooles.Solve() ;
#endif 
#ifdef TEST_KUNDERT
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
  }

  SparseDirectTimingVars::SS_Result.Set_Total_Time( TotalTime.ElapsedTime() ); 

  //
  //  Compute the error = norm(xcomp - xexact )
  //
  vector <double> error(numsolves) ; 
  double max_error = 0.0;
  
  passtemp->Update(1.0, *passx, -1.0, *passxexact, 0.0);

  passtemp->Norm2(&error[0]);
  for ( int i = 0 ; i< numsolves; i++ ) 
    if ( error[i] > max_error ) max_error = error[i] ; 
  SparseDirectTimingVars::SS_Result.Set_Error(max_error) ;

  //  passxexact->Norm2(&error[0] ) ; 
  //  passx->Norm2(&error ) ; 

  //
  //  Compute the residual = norm(Ax - b)
  //
  vector <double> residual(numsolves) ; 
  double max_resid = 0.0;
  
  passA->Multiply( transpose, *passx, *passtmp);
  passtemp->Update(1.0, *passtmp, -1.0, *passb, 0.0); 
  passtemp->Norm2(&residual[0]);

  for ( int i = 0 ; i< numsolves; i++ ) 
    if ( residual[i] > max_resid ) max_resid = residual[i] ; 


  SparseDirectTimingVars::SS_Result.Set_Residual(max_resid) ;
    
  vector <double> bnorm(numsolves); 
  passb->Norm2( &bnorm[0] ) ; 
  SparseDirectTimingVars::SS_Result.Set_Bnorm(bnorm[0]) ;

  vector <double> xnorm(numsolves); 
  passx->Norm2( &xnorm[0] ) ; 
  SparseDirectTimingVars::SS_Result.Set_Xnorm(xnorm[0]) ;

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;

  Comm.Barrier();
   return 0;
}
