#include "Epetra_CrsMatrix.h"
#include "TestAllClasses.h"
#include "TestOtherClasses.h"
#include "TestSuperludist.h"
 
int TestAllClasses(Epetra_CrsMatrix *& Amat, 
		   bool transpose, 
		   bool verbose, 
		   int Levels,
		   const double Rcond,
		   double &maxrelerror, 
		   double &maxrelresidual,
		   int &NumTests) {

  int errors = 0 ;

  if ( verbose) cout << " Testing Superludist " << endl ; 
  
  errors += TestSuperludist(Amat, 
			    transpose, 
			    verbose, 
			    Levels, 
			    Rcond, 
			    maxrelerror, 
			    maxrelresidual, 
			    NumTests ) ;
  
  if ( verbose) cout << " Testing UMFPACK " << endl ; 

  errors += TestOtherClasses(AMESOS_UMFPACK,
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing DSCPACK " << endl ; 

  errors += TestOtherClasses(AMESOS_DSCPACK,
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing KLU " << endl ; 

  errors += TestOtherClasses(AMESOS_KLU,
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing MUMPS " << endl ; 

  errors += TestOtherClasses(AMESOS_MUMPS,
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;


  if ( verbose) cout << " TestAllClasses errors = " << errors << endl ; 

  return errors;
}

