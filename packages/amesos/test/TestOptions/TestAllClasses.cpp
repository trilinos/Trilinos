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

  errors += TestOtherClasses("Amesos_Umfpack",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing DSCPACK " << endl ; 

  errors += TestOtherClasses("Amesos_Dscpack",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing SCALAPACK " << endl ; 

  errors += TestOtherClasses("Amesos_Scalapack",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing KLU " << endl ; 

  errors += TestOtherClasses("Amesos_Klu",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing MUMPS " << endl ; 

  errors += TestOtherClasses("Amesos_Mumps",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;

  if ( verbose) cout << " Testing SUPERLU " << endl ; 

  errors += TestOtherClasses("Amesos_Superlu",
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

