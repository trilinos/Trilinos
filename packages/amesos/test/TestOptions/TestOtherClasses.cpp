#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "TestOtherClasses.h"
#include "PerformOneSolveAndTest.h"

int TestOtherClasses(char* AmesosClass,
		     Epetra_CrsMatrix *& Amat, 
		     bool transpose, 
		     bool verbose, 
		     int Levels,
		     const double Rcond,
		     double &maxrelerror, 
		     double &maxrelresidual,
		     int &NumTests ) {

  int NumErrors = 0 ;
  maxrelerror = 0.0;
  maxrelresidual = 0.0;
  const Epetra_Comm& Comm = Amat->Comm();

  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", true );
    if ( verbose ) ParamList.set( "DebugLevel", 1 );

    //  ParamList.print( cerr, 10 ) ; 

    double relerror;
    double relresidual;
   
    int Errors = PerformOneSolveAndTest(AmesosClass,
					Comm, 
					transpose, 
					verbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual ) ; 

    if (Errors < 0 ) {
      if (verbose ) cout << AmesosClass << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 

      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestOtherClasses relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestOtherClasses relerror = " << relerror << endl ; 
      if (verbose) cout << " TestOtherClasses maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestOtherClasses maxrelerror = " << maxrelerror << endl ; 

    }
    if (verbose)  cout << " TestOtherClasses NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
    return NumErrors; 
  }
}
