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
#include "CreateTridi.h"
//
//  PartialFactorization checks to make sure that we clean up after our mess
//  before everything is done.
//  

const int MaxNumSteps = 7 ; 

int PartialFactorizationOneStep( char* AmesosClass,
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

  const Epetra_Map *Map = &Amat->RowMap() ; 

  Epetra_Vector xexact(*Map);
  Epetra_Vector x(*Map);

  Epetra_Vector b(*Map);
  Epetra_Vector bcheck(*Map);

  Epetra_Vector difference(*Map);

  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos Afactory;

  Abase = Afactory.Create( AmesosClass, Problem ) ; 

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
    if (verbose) ParamList.set( "DebugLevel", 1 );
    if (verbose) ParamList.set( "OutputLevel", 1 );
    if ( Steps > 1 ) {
      OUR_CHK_ERR( Abase->SetParameters( ParamList ) ); 
      if ( Steps > 2 ) {
		
	int ind[1];
	double val[1];
	ind[0] = 0;
	xexact.Random();
	xexact.PutScalar(1.0);
	
	//
	//  Compute cAx = A' xexact
	//
	double Value = 1.0 ;
	
	Amat->Multiply( transpose, xexact, b ) ;  //  b = A x2 = A A' A'' xexact
	
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
		  if (verbose) cout << " norm( x - xexact ) / norm(x) = " 
				    << norm_diff /norm_one << endl ; 
		  
		  relerror = norm_diff / norm_one ; 
		  
		  relresidual = norm_diff / norm_one ; 
		  
		  if (iam == 0 ) {
		    if ( relresidual * Rcond > 1e-16 ) {
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


int PartialFactorization( char* AmesosClass,
			  const Epetra_Comm &Comm, 
			  bool transpose, 
			  bool verbose, 
			  Teuchos::ParameterList ParamList, 
			  Epetra_CrsMatrix *& Amat, 
			  double Rcond ) {

#if 1
  for( int i =0 ; i < MaxNumSteps ; i ++ ) {
    string AC = AmesosClass ; 


    //  Amesos_Dscpack dies if SymbolicFactorization() is called twice in a row - Bug #1237 
    //  Hence we do not test Amesos_Dscpack with 3 or more steps here.
    if ( AC != "Amesos_Dscpack" || i < 3 ) { 
      PartialFactorizationOneStep( AmesosClass,
				   Comm, 
				   transpose, 
				   verbose, 
				   ParamList, 
				   Amat, 
				   Rcond,
				   i ) ;
    }
  }
#endif
}
