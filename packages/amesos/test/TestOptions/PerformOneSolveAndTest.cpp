//
//  OUR_CHK_ERR always returns 1 on error.
//

#define OUR_CHK_ERR(a) { { int epetra_err = a; \
                      if (epetra_err != 0) { cerr << "Amesos ERROR " << epetra_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; \
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
			    double& relresidual )
{

  RefCountPtr<Epetra_CrsMatrix> MyMat ; 
  RefCountPtr<Epetra_CrsMatrix> MyMatWithDiag ; 

  PartialFactorization( AmesosClass, Comm, transpose, verbose, ParamList, InMat, Rcond );

  MyMat = rcp( new Epetra_CrsMatrix( *InMat ) ); 

  Amesos_TestRowMatrix ATRW( &*MyMat ) ; 
  
  Epetra_RowMatrix* MyRowMat ; 
  
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
    break;
  }

  if (ParamList.isParameter("AddToDiag")) { 
    //
    //  If AddToDiag is set, create a matrix which is numerically identical, but structurally 
    //  has no missing diagaonal entries.   In other words, every diagonal element in MyMayWithDiag 
    //  has an entry in the matrix, though that entry will be zero if InMat has no entry for that
    //  particular diagonal element.  
    //
    MyMatWithDiag = NewMatNewMap( *InMat, 2, 0, 0, 0, 0 );  //  Ensure that all diagonal entries exist ;

    //
    //  Now add AddToDiag to each diagonal element.  
    //
    double AddToDiag = ParamList.get("AddToDiag", AddToDiag );
    Epetra_Vector Diag( MyMatWithDiag->RowMap() );
    Epetra_Vector AddConstVecToDiag( MyMatWithDiag->RowMap() );
    AddConstVecToDiag.PutScalar( AddToDiag );
    assert( MyMatWithDiag->ExtractDiagonalCopy( Diag ) == 0 );
    Diag.Update( 1.0, AddConstVecToDiag, 1.0 ) ; 
    assert(MyMatWithDiag->ReplaceDiagonalValues( Diag ) == 0 ) ; 
  } else { 
    MyMatWithDiag = rcp( new Epetra_CrsMatrix( *InMat ) ); 
  }
  //  Epetra_CrsMatrix*& Amat = &*MyMat ; 

  if ( verbose ) cout << " Partial Factorization complete " << endl ; 

  relerror = 0 ; 
  relresidual = 0 ; 

  assert( Levels >= 1 && Levels <= 3 ) ; 

  int iam = Comm.MyPID() ; 
  int errors = 0 ; 

  const Epetra_Map *Map = &MyMat->RowMap() ; 
  const Epetra_Map *RangeMap = &MyMat->OperatorRangeMap() ; 
  const Epetra_Map *DomainMap = &MyMat->OperatorDomainMap() ; 

  Epetra_Vector xexact(*DomainMap);
  Epetra_Vector x(*DomainMap);

  Epetra_Vector cAx(*DomainMap);
  Epetra_Vector sAx(*DomainMap);
  Epetra_Vector kAx(*DomainMap);

  Epetra_Vector cAAx(*DomainMap);
  Epetra_Vector sAAx(*DomainMap);
  Epetra_Vector kAAx(*DomainMap);

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
    Epetra_CrsMatrix* ECM = dynamic_cast<Epetra_CrsMatrix*>(MyRowMat) ; 
    
    //
    //  We only set transpose if we have to - this allows valgrind to check
    //  that transpose is set to a default value before it is used.
    //
    if ( transpose ) OUR_CHK_ERR( Abase->SetUseTranspose( transpose ) ); 
    if (verbose) ParamList.set( "DebugLevel", 1 );
    if (verbose) ParamList.set( "OutputLevel", 1 );
    OUR_CHK_ERR( Abase->SetParameters( ParamList ) ); 
#if 0 
    OUR_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
    OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 
#endif 

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

    if ( verbose ) cout << " Compute  b = A x2 = A A' A'' xexact  " << endl ; 

    MyMatWithDiag->Multiply( transpose, cAAx, b ) ;  //  b = A x2 = A A' A'' xexact
 

    //  Phase 2:  Solve A' A' A x = b 
    //
    //
    //  Solve A sAAx = b 
    //
    Problem.SetLHS( &sAAx );
    Problem.SetRHS( &b );


    OUR_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
    OUR_CHK_ERR( Abase->SymbolicFactorization(  ) );     // This should be irrelevant, but should nonetheless be legal 
    OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 
    OUR_CHK_ERR( Abase->Solve(  ) ); 

    if ( Levels >= 2 ) 
      {
        OUR_CHK_ERR( Abase->SetUseTranspose( transpose ) ); 
	Problem.SetLHS( &sAx );
	Problem.SetRHS( &sAAx );
	val[0] =  Value ; 
	if ( MyMat->MyGRID( 0 ) )
	  MyMat->SumIntoMyValues( 0, 1, val, ind ) ; 
	OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 
	
	Teuchos::ParameterList* NullList = (Teuchos::ParameterList*) 0 ;  
	//      We do not presently handle null lists.
	//	OUR_CHK_ERR( Abase->SetParameters( *NullList ) );   // Make sure we handle null lists 
	OUR_CHK_ERR( Abase->Solve(  ) ); 
	
      }
    else
      {
	sAx = sAAx ;
      }

    if ( Levels >= 3 ) 
      {
	Problem.SetLHS( &x );
	Problem.SetRHS( &sAx );
	OUR_CHK_ERR( Abase->Solve(  ) ); 
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
	    cout << " TestOptions requires a non-zero entry in A(1,1) " << endl ; 
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

    if ( verbose ) cout << " Levels =  " << Levels << endl ; 
    if ( verbose ) cout << " Rcond =  " << Rcond << endl ; 

    double norm_diff ;
    double norm_one ;

    DomainDiff.Update( 1.0, sAAx, -1.0, cAAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAAx.Norm2( &norm_one ) ; 
    if (verbose) cout << " norm( sAAx - cAAx ) / norm(sAAx ) = " 
		      << norm_diff /norm_one << endl ; 


    DomainDiff.Update( 1.0, sAx, -1.0, cAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAx.Norm2( &norm_one ) ; 
    if (verbose) cout << " norm( sAx - cAx ) / norm(sAx ) = " 
		      << norm_diff /norm_one << endl ; 


    DomainDiff.Update( 1.0, x, -1.0, xexact, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    x.Norm2( &norm_one ) ; 
    if (verbose) cout << " norm( x - xexact ) / norm(x) = " 
		      << norm_diff /norm_one << endl ; 

    relerror = norm_diff / norm_one ; 

    DomainDiff.Update( 1.0, sAx, -1.0, kAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAx.Norm2( &norm_one ) ; 
    if (verbose) cout << " norm( sAx - kAx ) / norm(sAx ) = " 
		      << norm_diff /norm_one << endl ; 


    DomainDiff.Update( 1.0, sAAx, -1.0, kAAx, 0.0 ) ;
    DomainDiff.Norm2( &norm_diff ) ; 
    sAAx.Norm2( &norm_one ) ; 
    if (verbose) cout << " norm( sAAx - kAAx ) / norm(sAAx ) = " 
		      << norm_diff /norm_one << endl ; 


    RangeDiff.Update( 1.0, bcheck, -1.0, b, 0.0 ) ;
    RangeDiff.Norm2( &norm_diff ) ; 
    bcheck.Norm2( &norm_one ) ; 
    if (verbose) cout << " norm( bcheck - b ) / norm(bcheck ) = " 
		      << norm_diff /norm_one << endl ; 

    relresidual = norm_diff / norm_one ; 


    if (iam == 0 ) {
      if ( relresidual * Rcond < 1e-16 ) {
	if (verbose) cout << " Test 1 Passed " << endl ;
      } else {
	if (verbose) cout << " TEST 1 FAILED " << endl ;
	errors += 1 ; 
      }
    }


    //  #define FACTOR_B
#ifdef FACTOR_B


    garbage;

    //
    //  Now we check to make sure that we can change the problem and 
    //  re-factorize.  
    //




    int BNumPoints;
    const Epetra_Map* BMap ; 
    string Aclass = AmesosClass ;
    const bool AllowDiffProbSize = ( Aclass != "Amesos_Superludist" ) ;
    if ( AllowDiffProbSize ) { 
      BNumPoints = 6;  // Must be between 2 and 100 (on large matrices,
      // the problem is quite ill-conditioned) 
    
      // Construct a Map that puts approximately the same number of 
      // equations on each processor.
      BMap = new Epetra_Map(BNumPoints, 0, Comm);
    }  else {
      BMap = Map ;
      BNumPoints = BMap->NumGlobalElements();
    }
    
      
    //  Create an empty EpetraCrsMatrix 
    Epetra_CrsMatrix B(Copy, (*BMap), 0);

    //
    //  Populate A with a [-1,2,-1] tridiagonal matrix WITH -1 in the
    //  off diagonal corners.
    //  See CreateTridi.cpp in this directory 
    CreateTridiPlus( B ) ; 

    Epetra_Vector Bx2((*BMap)), Bx1((*BMap)), Bx((*BMap)), Bb((*BMap)), Bresidual((*BMap)), Btemp((*BMap));

    //

    //
    //  Factor B
    //
    Problem.SetOperator( &B );
    OUR_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
    OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 

    Bb.Random();
    Bb.PutScalar(1.0);
    //
    //  Solve B x = b 
    //
    Problem.SetLHS( &Bx );
    Problem.SetRHS( &Bb );
    OUR_CHK_ERR( Abase->Solve(  ) ); 


    ind[0] = 0;
    val[0] = 1 ; 
    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 

    //  if (verbose) cout << " B' = " << B << endl ; 
    OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 

    //
    //  Solve B' x1 = x 
    //
    Problem.SetLHS( &Bx1 );
    Problem.SetRHS( &Bx );
    OUR_CHK_ERR( Abase->Solve(  ) ); 

    //  if (verbose) cout << " x1 = " << x1 << endl ; 

    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 

    //  if (verbose) cout << " B'' = " << B << endl ; 
    OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 

    //
    //  Solve B" x2 = x1
    //
    Problem.SetLHS( &Bx2 );
    Problem.SetRHS( &Bx1 );
    OUR_CHK_ERR( Abase->Solve(  ) ); 

    //  if (verbose) cout << " x2 = " << x2 << endl ; 

    //
    //  Compute the residual: B B' B" x2 - b
    //

    B.Multiply( transpose, Bx2, Btemp ) ; //  temp = B x2

    //  if (verbose) cout << " temp = " << temp << endl ; 

    val[0] = -val[0] ; 
    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 
    B.Multiply( transpose, Btemp, Bx2 ) ; //  x2 = B' B" x2



    //  if (verbose) cout << " x2 = " << x2 << endl ; 


    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 
    B.Multiply( transpose, Bx2, Btemp ) ; //  temp = B B' B'' x2


    //  if (verbose) cout << " temp = " << temp << endl ; 
    //  if (verbose) cout << " b = " << b << endl ; 




    Bresidual.Update( 1.0, Btemp, -1.0, Bb, 0.0 ) ;
    //  if (verbose) cout << " residual = " << residual << endl ; 

    double norm_residual;
    Bresidual.Norm2( &norm_residual ) ; 

    if (iam == 0 ) {
      if (verbose) cout << " norm2(B B' B'' x-b) = " << norm_residual << endl ; 
      //
      //  This is an ill-conditioned problem
      //
      if ( norm_residual < (1e-15)*(1.0*BNumPoints*BNumPoints*BNumPoints*
				    BNumPoints*BNumPoints*BNumPoints) ) {
	if (verbose) cout << " Test 2 Passed " << endl ;
      } else {
	if (verbose) cout << " TEST 2 FAILED " << endl ;
	errors += 1 ; 
      }
    }
    if ( AllowDiffProbSize ) {
      delete BMap;
    }
#endif

    delete Abase;
  }

  return errors;
  
}


