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
//
//  PartialFactorization checks to make sure that we clean up after our mess
//  before everything is done.
//  

const int MaxNumSteps = 7 ; 

int PartialFactorizationOneStep( const char* AmesosClass,
				 const Epetra_Comm &Comm, 
				 bool transpose, 
				 bool verbose, 
				 Teuchos::ParameterList ParamList, 
				 Epetra_CrsMatrix *& Amat, 
				 double Rcond, 
				 int Steps ) 
{
	
  assert( Steps >= 0 && Steps < MaxNumSteps ) ; 

  int iam = Comm.MyPID() ; 
  int errors = 0 ; 

  const Epetra_Map *RangeMap = 
    transpose?&Amat->OperatorDomainMap():&Amat->OperatorRangeMap() ; 
  const Epetra_Map *DomainMap = 
    transpose?&Amat->OperatorRangeMap():&Amat->OperatorDomainMap() ; 

  Epetra_Vector xexact(*DomainMap);
  Epetra_Vector x(*DomainMap);

  Epetra_Vector b(*RangeMap);
  Epetra_Vector bcheck(*RangeMap);

  Epetra_Vector difference(*DomainMap);

  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos Afactory;

  Abase = Afactory.Create( AmesosClass, Problem ) ; 

  std::string AC = AmesosClass ;
  if ( AC == "Amesos_Mumps" ) { 
    ParamList.set( "NoDestroy", true );
   Abase->SetParameters( ParamList ) ; 
  }

  double relerror = 0 ; 
  double relresidual = 0 ; 
  
  if ( Steps > 0 ) {
    //
    //  Phase 1:  Compute b = A' A' A xexact
    //
    Problem.SetOperator( Amat );
   
    //
    //  We only set transpose if we have to - this allows valgrind to check
    //  that transpose is set to a default value before it is used.
    //
    if ( transpose ) OUR_CHK_ERR( Abase->SetUseTranspose( transpose ) ); 
    //    if (verbose) ParamList.set( "DebugLevel", 1 );
    //    if (verbose) ParamList.set( "OutputLevel", 1 );
    if ( Steps > 1 ) {
      OUR_CHK_ERR( Abase->SetParameters( ParamList ) ); 
      if ( Steps > 2 ) {
		
	int ind[1];
	ind[0] = 0;
	xexact.Random();
	xexact.PutScalar(1.0);
	
	//
	//  Compute cAx = A' xexact
	//
	Amat->Multiply( transpose, xexact, b ) ;  //  b = A x2 = A A' A'' xexact

#if 0 
	std::cout << __FILE__ << "::"  << __LINE__ << "b = " << std::endl ; 
	b.Print( std::cout ) ; 
	std::cout << __FILE__ << "::"  << __LINE__ << "xexact = " << std::endl ; 
	xexact.Print( std::cout ) ; 
	std::cout << __FILE__ << "::"  << __LINE__ << "x = " << std::endl ; 
	x.Print( std::cout ) ; 
#endif
	//
	//  Phase 2:  Solve A' A' A x = b 
	//
	//
	//  Solve A sAAx = b 
	//
	Problem.SetLHS( &x );
	Problem.SetRHS( &b );
	OUR_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
	if ( Steps > 2 ) {
	  OUR_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
	  if ( Steps > 3 ) {
	    OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 
	    if ( Steps > 4 ) {
	      OUR_CHK_ERR( Abase->NumericFactorization(  ) ); 
	      if ( Steps > 5 ) {
		OUR_CHK_ERR( Abase->Solve(  ) ); 
		if ( Steps > 6 ) {
		  OUR_CHK_ERR( Abase->Solve(  ) ); 


		  Amat->Multiply( transpose, x, bcheck ) ; //  temp = A" x2
		  
		  double norm_diff ;
		  double norm_one ;
		  
		  difference.Update( 1.0, x, -1.0, xexact, 0.0 ) ;
		  difference.Norm2( &norm_diff ) ; 
		  x.Norm2( &norm_one ) ; 
		  
		  relerror = norm_diff / norm_one ; 
		  
		  relresidual = norm_diff / norm_one ; 
		  
		  if (iam == 0 ) {
		    if ( relresidual * Rcond > 1e-16 ) {
		      if (verbose) std::cout << __FILE__ << "::"<< __LINE__ 
					<< " norm( x - xexact ) / norm(x) = " 
					<< norm_diff /norm_one << std::endl ; 
		      errors += 1 ; 
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
}
 delete Abase;
 
 return errors;
 
}


int PartialFactorization( const char* AmesosClass,
			  const Epetra_Comm &Comm, 
			  bool transpose, 
			  bool verbose, 
			  Teuchos::ParameterList ParamList, 
			  Epetra_CrsMatrix *& Amat, 
			  double Rcond ) {

  double relerror = 0 ; 
  double relresidual = 0 ; 

  for( int i =0 ; i < MaxNumSteps ; i ++ ) {
    std::string AC = AmesosClass ; 

    //
    //  I have turned this test back on because it no longer seems to fail.  Although Amesos_Dscpack
    //  fails other aspects of TestOptions.
    //
    //  Amesos_Dscpack dies if SymbolicFactorization() is called twice in a row - Bug #1237 
    //  Hence we do not test Amesos_Dscpack with 3 or more steps here.
    //    if ( AC != "Amesos_Dscpack" || i < 3 ) { 
    {
      OUR_CHK_ERR( PartialFactorizationOneStep( AmesosClass,
						   Comm, 
						   transpose, 
						   verbose, 
						   ParamList, 
						   Amat, 
						   Rcond,
						   i ) ) ;
    }
  }
  return(0);
}
