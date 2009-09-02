//
//  OUR_CHK_ERR always returns 1 on error.
//

#define OUR_CHK_ERR(a) { { int epetra_err = a; \
                      if (epetra_err != 0) { std::cerr << "Amesos ERROR " << epetra_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << std::endl; \
relerror = 1.3e15; relresidual=1e15; return(1);}  }\
                   }

#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "PerformOneSolveAndTest.h"
#include "PartialFactorization.h"
#include "CreateTridi.h"
#include "NewMatNewMap.h" 
#include "Amesos_TestRowMatrix.h" 

//
//  Returns the number of failures.
//  Note:  If AmesosClass is not supported, PerformOneSolveAndTest() will 
//  always return 0
//
//  Still have to decide where we are going to check the residual.  
//
//  The following table shows the variable names that we use for 
//  each of the three phases:  
//     compute - which computes the correct value of b
//     solve - which solves for x in  A' A' A x = b 
//     check - which computes bcheck = A' A' A x 
//
//  For ill-conditioned matrices we restrict the test to one or two 
//  solves, by setting Levels to 1 or 2 on input to this routine.
//  When Levels is less than 3, some of the transformations
//  shown in the table as "->" and "<-" are not performed, instead 
//  a direct copy is made.
//
//  In the absence of roundoff, each item in a given column should 
//  be identical.  
//
//  If Levels = 3, we compute and solve A' A' A x = b and hence all 
//  of the transformations are performed
//
//  If Levels = 2, the transformations shown in the first column of 
//  transformations (labelled Levels>=3) are replaced by a direct copy.
//
//  If Levels = 1, only the transformations shown in the third column
//  are performed, the others being replaced by direct copies.
//  
//                           Levels>=3    Levels>=2
//                              A'         A'            A
//  compute             xexact  ->  cAx    ->     cAAx   ->       b 
//  solve               x       <-  sAx    <-     sAAx   <-       b
//  check               x       ->  kAx    ->     kAAx   ->  bcheck
//
//  Note that since Levels 2 and 3 use the same A, there is no need to 
//  call NumericFactorization() between the second and third call to Solve. 
//   
//  On testing AddToDiag
//  ====================
//
//  When this code was written, our intent was to have AddToDiag add a constant 
//  value to the diagonal even for non-existent diagonal elements.  For now, we have
//  backed off of that.  If we decide to go back to the earlier semantics for 
//  AddToDiag (i.e. that it would force all diagonal elements to exist in the 
//  matrix, this can be tested by setting AddToAllDiagonalElements to true ; 
//
//  See bug #1928 for further discussion.
//
//  ExpectedError
//  =============
//
//  ExpectedError allows us to test for Singular and Structurally Singular matrices
//    ExpectedError = StructurallySingularMatrix or NumericallySingularMatrix


int PerformOneSolveAndTest( const char* AmesosClass,
			    int EpetraMatrixType,
			    const Epetra_Comm &Comm, 
			    bool transpose, 
			    bool verbose, 
			    Teuchos::ParameterList ParamList, 
			    Epetra_CrsMatrix *& InMat, 
			    int Levels, 
			    const double Rcond,
			    double& relerror,
			    double& relresidual,
			    int ExpectedError )
{

  bool AddToAllDiagonalElements =  ParamList.get( "AddZeroToDiag", false ) ;


  bool TrustMe = ParamList.get( "TrustMe", false );

  bool MyVerbose = false ; //  setting this equal to verbose produces too much output and exceeds the test harnesses 1 Megabyte limit

  RCP<Epetra_CrsMatrix> MyMat ; 
  RCP<Epetra_CrsMatrix> MyMatWithDiag ; 

  MyMat = rcp( new Epetra_CrsMatrix( *InMat ) ); 

  Amesos_TestRowMatrix ATRW( &*MyMat ) ; 
  
  Epetra_RowMatrix* MyRowMat = 0; 
  
  assert ( EpetraMatrixType >= 0 && EpetraMatrixType <= 2 );
  switch ( EpetraMatrixType ) {
  case 0:
    MyRowMat = &*MyMat ; 
    break;
  case 1:
    MyRowMat = &ATRW;
    break;
  case 2:
    MyMat->OptimizeStorage(); 
    MyRowMat = &*MyMat ; 
    bool OptStorage = MyMat->StorageOptimized();
    assert( OptStorage) ; 
    break;
  }
  bool OptStorage = MyMat->StorageOptimized();
  
  Epetra_CrsMatrix* MatPtr = &*MyMat ;

  const std::string AC = AmesosClass ;
  if ( ExpectedError == 0 ) { 
    if ( AC != "Amesos_Pardiso" ) {
      OUR_CHK_ERR ( PartialFactorization( AmesosClass, Comm, transpose, MyVerbose, 
					  ParamList, MatPtr, Rcond ) );
    } else {
      if (MyVerbose) std::cout << " AC = "  << AC << " not tested in Partial Factorization" <<std::endl ;   // bug #1915
    }
  }

  if (ParamList.isParameter("AddToDiag") ) { 
      int oldtracebackmode = InMat->GetTracebackMode( ) ;   

    //
    //  If AddToDiag is set, create a matrix which is numerically identical, but structurally 
    //  has no missing diagaonal entries.   In other words, every diagonal element in MyMayWithDiag 
    //  has an entry in the matrix, though that entry will be zero if InMat has no entry for that
    //  particular diagonal element.  
    //
    //  It turns out that this does not actually make sure that all diagonal entries exist because it 
    //  does not deal with maps that are missing the diagonal element.
    //
    if ( AddToAllDiagonalElements ) { 
      MyMatWithDiag = NewMatNewMap( *InMat, 2, 0, 0, 0, 0 );  //  Ensure that all diagonal entries exist ;
    } else { 
      InMat->SetTracebackMode( EPETRA_MIN(oldtracebackmode,1) ) ;   // If the matrix is mimssing diagonal elements, the call to ReplaceDiagonalElements will return 1 indicating missing diagonal elements
      MyMatWithDiag = NewMatNewMap( *InMat, 0, 0, 0, 0, 0 );  //  Leave the matrix unchanged structurally 
    }
    //
    //  Now add AddToDiag to each diagonal element.  
    //
    double AddToDiag = ParamList.get("AddToDiag", 0.0 );
    Epetra_Vector Diag( MyMatWithDiag->RowMap() );
    Epetra_Vector AddConstVecToDiag( MyMatWithDiag->RowMap() );
    AddConstVecToDiag.PutScalar( AddToDiag );

    assert( MyMatWithDiag->ExtractDiagonalCopy( Diag ) == 0 );
    Diag.Update( 1.0, AddConstVecToDiag, 1.0 ) ; 
    assert(MyMatWithDiag->ReplaceDiagonalValues( Diag ) >= 0 ) ;   // This may return 1 indicating that structurally non-zero elements were left untouched. 

      InMat->SetTracebackMode( oldtracebackmode ) ;   
  } else { 
    MyMatWithDiag = rcp( new Epetra_CrsMatrix( *InMat ) ); 
  }

  if ( MyVerbose ) std::cout << " Partial Factorization complete " << std::endl ; 

  relerror = 0 ; 
  relresidual = 0 ; 

  assert( Levels >= 1 && Levels <= 3 ) ; 

  int iam = Comm.MyPID() ; 
  int errors = 0 ; 

  const Epetra_Map *RangeMap = 
    transpose?&MyMat->OperatorDomainMap():&MyMat->OperatorRangeMap() ; 
  const Epetra_Map *DomainMap = 
    transpose?&MyMat->OperatorRangeMap():&MyMat->OperatorDomainMap() ; 

  Epetra_Vector xexact(*DomainMap);
  Epetra_Vector x(*DomainMap);

  Epetra_Vector cAx(*DomainMap);
  Epetra_Vector sAx(*DomainMap);
  Epetra_Vector kAx(*DomainMap);

  Epetra_Vector cAAx(*DomainMap);
  Epetra_Vector sAAx(*DomainMap);
  Epetra_Vector kAAx(*DomainMap);

  Epetra_Vector FixedLHS(*DomainMap);
  Epetra_Vector FixedRHS(*RangeMap);

  Epetra_Vector b(*RangeMap);
  Epetra_Vector bcheck(*RangeMap);

  Epetra_Vector DomainDiff(*DomainMap);
  Epetra_Vector RangeDiff(*RangeMap);

  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos Afactory;




  Abase = Afactory.Create( AmesosClass, Problem ) ; 

  relerror = 0 ; 
  relresidual = 0 ; 

  if ( Abase == 0 ) 
    assert( false ) ; 
  else {

    //
    //  Phase 1:  Compute b = A' A' A xexact
    //
    Problem.SetOperator( MyRowMat );
    //    Epetra_CrsMatrix* ECM = dynamic_cast<Epetra_CrsMatrix*>(MyRowMat) ; 

    //
    //  We only set transpose if we have to - this allows valgrind to check
    //  that transpose is set to a default value before it is used.
    //
    if ( transpose ) OUR_CHK_ERR( Abase->SetUseTranspose( transpose ) ); 
    if (MyVerbose) ParamList.set( "DebugLevel", 1 );
    if (MyVerbose) ParamList.set( "OutputLevel", 1 );
    OUR_CHK_ERR( Abase->SetParameters( ParamList ) ); 

    if ( TrustMe ) {
      Problem.SetLHS( &FixedLHS );
      Problem.SetRHS( &FixedRHS );
      assert( OptStorage) ;
    }

    // bug #2184
    //  Structurally singular matrices are not detected by 
    //  Amesos_Klu::SymbolicFactorization() but they are detected by
    //  Amesos_Klu::NumericFactorization()
    if ( ExpectedError == StructurallySingularMatrixError )
      ExpectedError = NumericallySingularMatrixError ;
    if ( ExpectedError == StructurallySingularMatrixError ) {
      Epetra_CrsMatrix* ECM = dynamic_cast<Epetra_CrsMatrix*>(MyRowMat) ; 
      int oldtracebackmode = ECM->GetTracebackMode( ) ;         
      ECM->SetTracebackMode(0);  // We expect an error, but we don't want it to print out

      const int SymbolicFactorizationReturn = Abase->SymbolicFactorization(  ) ; 
      ECM->SetTracebackMode(oldtracebackmode);
      // bug #2245 - Amesos fails to return error consistently across all 
      // processes.  When this bug is fixed, remove "iam == 0 &&" from the next line 
      if ( iam == 0 && SymbolicFactorizationReturn != ExpectedError ) {
	std::cout << " SymbolicFactorization returned " << SymbolicFactorizationReturn 
	     << " should be " << ExpectedError << std::endl ; 
	OUR_CHK_ERR( 1 ) ; 
      } else { 
	return 0;   //  Returned the correct error code for this matrix 
      }
    } else { 
      const int SymbolicFactorizationReturn = Abase->SymbolicFactorization(  ) ; 
      OUR_CHK_ERR( SymbolicFactorizationReturn ) ; 
    }
    if ( ExpectedError == NumericallySingularMatrixError ) {
      Epetra_CrsMatrix* ECM = dynamic_cast<Epetra_CrsMatrix*>(MyRowMat) ; 
      int oldtracebackmode = ECM->GetTracebackMode( ) ;   
      ECM->SetTracebackMode(0); // We expect an error, but we don't want it to print out

      const int NumericFactorizationReturn = Abase->NumericFactorization(  ) ; 
      ECM->SetTracebackMode(oldtracebackmode);
      // bug #2245 - Amesos fails to return error consistently across all 
      // processes.  When this bug is fixed, remove "iam == 0 &&" from the next line 
      if ( iam == 0 && NumericFactorizationReturn != ExpectedError ) {
	std::cout << " NumericFactorization returned " << NumericFactorizationReturn 
	     << " should be " << ExpectedError << std::endl ; 
	OUR_CHK_ERR( 1 ) ; 
      } else { 
	return 0;   //  Returned the correct error code for this matrix 
      }
    } else {
      const int NumericFactorizationReturn = Abase->NumericFactorization(  ) ; 
      OUR_CHK_ERR( NumericFactorizationReturn ) ; 
    }
    int ind[1];
    double val[1];
    ind[0] = 0;
    xexact.Random();
    xexact.PutScalar(1.0);

    //
    //  Compute cAx = A' xexact
    //
    double Value = 1.0 ;
    if ( Levels == 3 ) 
      {
	val[0] = Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) ) { 
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
	}
	MyMatWithDiag->Multiply( transpose, xexact, cAx ) ; 

	val[0] = - Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) )
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
      }
    else
      {
	cAx = xexact ;
      }

    //
    //  Compute cAAx = A' cAx
    //
    if ( Levels >= 2 ) 
      {
	val[0] =  Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) )
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
	MyMatWithDiag->Multiply( transpose, cAx, cAAx ) ; //  x2 = A' x1

	val[0] = - Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) )
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
      }
    else
      {
	cAAx = cAx ;
      }

    if ( MyVerbose ) std::cout << " Compute  b = A x2 = A A' A'' xexact  " << std::endl ; 

    MyMatWithDiag->Multiply( transpose, cAAx, b ) ;  //  b = A x2 = A A' A'' xexact
 

    //  Phase 2:  Solve A' A' A x = b 
    //
    //
    //  Solve A sAAx = b 
    //
    if ( TrustMe ) { 
      FixedRHS = b;
    } else { 
      Problem.SetLHS( &sAAx );
      Problem.SetRHS( &b );
    }


    OUR_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
    OUR_CHK_ERR( Abase->SymbolicFactorization(  ) );     // This should be irrelevant, but should nonetheless be legal 
    OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 
    OUR_CHK_ERR( Abase->Solve(  ) ); 
    if ( TrustMe ) sAAx = FixedLHS ; 

    if ( Levels >= 2 ) 
      {
        OUR_CHK_ERR( Abase->SetUseTranspose( transpose ) ); 
	if ( TrustMe ) { 
	  FixedRHS = sAAx ;
	} else { 
	  Problem.SetLHS( &sAx );
	  Problem.SetRHS( &sAAx );
	}
	val[0] =  Value ; 
	if ( MyMat->MyGRID( 0 ) )
	  MyMat->SumIntoMyValues( 0, 1, val, ind ) ; 
	OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 
	OUR_CHK_ERR( Abase->Solve(  ) ); 
	if ( TrustMe ) sAx = FixedLHS ; 
	
      }
    else
      {
	sAx = sAAx ;
      }

    if ( Levels >= 3 ) 
      {
	if ( TrustMe ) { 
	  FixedRHS = sAx ;
	} else { 
	  Problem.SetLHS( &x );
	  Problem.SetRHS( &sAx );
	}
	OUR_CHK_ERR( Abase->Solve(  ) ); 
	if ( TrustMe ) x = FixedLHS ;
      }
    else
      {
	x = sAx ; 
      }

    if ( Levels >= 2 ) 
      {
	val[0] =  -Value ; 
	if ( MyMat->MyGRID( 0 ) ) {
	  if ( MyMat->SumIntoMyValues( 0, 1, val, ind ) ) { 
	    std::cout << " TestOptions requires a non-zero entry in A(1,1) " << std::endl ; 
	  }
	}
      }

    //
    //  Phase 3:  Check the residual: bcheck = A' A' A x 
    //


    if ( Levels >= 3 ) 
      {
	val[0] =  Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) )
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
	MyMatWithDiag->Multiply( transpose, x, kAx ) ;
	val[0] =  -Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) )
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
      }
    else
      {
	kAx = x ; 
      }

    if ( Levels >= 2 ) 
      {
	val[0] =  Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) )
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
	MyMatWithDiag->Multiply( transpose, kAx, kAAx ) ;
	val[0] =  -Value ; 
	if ( MyMatWithDiag->MyGRID( 0 ) )
	  MyMatWithDiag->SumIntoMyValues( 0, 1, val, ind ) ; 
      }
    else
      {
	kAAx = kAx ; 
      }

    MyMatWithDiag->Multiply( transpose, kAAx, bcheck ) ; //  temp = A" x2

    double norm_diff ;
    double norm_one ;

    DomainDiff.Update( 1.0, sAAx, -1.0, cAAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAAx.Norm2( &norm_one ) ; 
    if (MyVerbose) std::cout << __FILE__ << "::" << __LINE__ 
			<< " norm( sAAx - cAAx ) / norm(sAAx ) = " 
			<< norm_diff /norm_one << std::endl ; 


    

    DomainDiff.Update( 1.0, sAx, -1.0, cAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAx.Norm2( &norm_one ) ; 
    if (MyVerbose) std::cout 
      << __FILE__ << "::" << __LINE__ 
      << " norm( sAx - cAx ) / norm(sAx ) = " 
		      << norm_diff /norm_one << std::endl ; 


    DomainDiff.Update( 1.0, x, -1.0, xexact, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    x.Norm2( &norm_one ) ; 
    if (MyVerbose) std::cout 
      << __FILE__ << "::" << __LINE__ 
      << " norm( x - xexact ) / norm(x) = " 
      << norm_diff /norm_one << std::endl ; 

    relerror = norm_diff / norm_one ; 

    DomainDiff.Update( 1.0, sAx, -1.0, kAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAx.Norm2( &norm_one ) ; 
    if (MyVerbose) std::cout 
      << __FILE__ << "::" << __LINE__ 
      << " norm( sAx - kAx ) / norm(sAx ) = " 
      << norm_diff /norm_one << std::endl ; 


    DomainDiff.Update( 1.0, sAAx, -1.0, kAAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAAx.Norm2( &norm_one ) ; 
    if (MyVerbose) std::cout 
      << __FILE__ << "::" << __LINE__ 
      << " norm( sAAx - kAAx ) / norm(sAAx ) = " 
      << norm_diff /norm_one << std::endl ; 


    RangeDiff.Update( 1.0, bcheck, -1.0, b, 0.0 ) ;
    RangeDiff.Norm2( &norm_diff ) ; 
    bcheck.Norm2( &norm_one ) ; 
    if (MyVerbose) std::cout 
      << __FILE__ << "::" << __LINE__ 
      << " norm( bcheck - b ) / norm(bcheck ) = " 
      << norm_diff /norm_one << std::endl ; 

    relresidual = norm_diff / norm_one ; 

    if (iam == 0 ) {
      if ( relresidual * Rcond < 1e-16 ) {
	if (MyVerbose) std::cout << " Test 1 Passed " << std::endl ;
      } else {
      std::cout <<  __FILE__ << "::"  << __LINE__ << 
	  " relresidual = " << relresidual <<
	  " TEST FAILED " <<
	  " ParamList = " << ParamList << std::endl ; 
	errors += 1 ; 
      }
    }

    delete Abase;
  }

  return errors;

}


