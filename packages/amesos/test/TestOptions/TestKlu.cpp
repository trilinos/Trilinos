#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "PerformOneSolveAndTest.h"
#include "TestKlu.h"

//
//  Tests:
//     1)  no parameters
//     2)  Refactorize = true 
//     3)  ScaleMethod = 1 - argh I don't see how ScaleMEthod can work
//     4)  ComputeTrueResidual==true
//
int TestKlu( Epetra_CrsMatrix *& Amat, 
	     int EpetraMatrixType,
	     const bool transpose, 
	     const bool verbose, 
	     const int Levels,
	     const double Rcond,
	     bool RowMapEqualsColMap, 
	     double &maxrelerror, 
	     double &maxrelresidual,
	     int &NumTests ) {
  
  int NumErrors = 0 ;
  maxrelerror = 0.0;
  maxrelresidual = 0.0;
  const Epetra_Comm& Comm = Amat->Comm();

  //
  //     1)  no parameters

  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 1) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
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
      NumErrors++;
      NumTests++ ; 
      if ( verbose ) {
	cout << "Amesos_Klu failed with error code " << Errors<< endl ; 
      }
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestKlu relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestKlu relerror = " << relerror << endl ; 
      if (verbose) cout << " TestKlu maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestKlu maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestKlu NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 1) no fail second time through yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
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
      NumErrors++;
      NumTests++ ; 
      if ( verbose ) {
	cout << "Amesos_Klu failed with error code " << Errors<< endl ; 
      }
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestKlu 2nd relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestKlu relerror = " << relerror << endl ; 
      if (verbose) cout << " TestKlu maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestKlu 2nd maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestKlu NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  //
  //     2)  Refactorize = true 
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
    ParamList.set( "Refactorize", true );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Klu Test 2) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual ) ; 
      
    if ( verbose ) cout << " Test 2) no fail here either " << endl ; 

    if (Errors < 0 ) {
      if (verbose ) cout << "Amesos_Klu" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestKlu relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestKlu relerror = " << relerror << endl ; 
      if (verbose) cout << " TestKlu maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestKlu maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestKlu NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ; 
    }
  }
  //
  //     2A)  AddToDiag 
  if (RowMapEqualsColMap ) {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
    ParamList.set( "Refactorize", true );
    ParamList.set( "AddToDiag", 1e-2 );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Klu Test 2) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual ) ; 
      
    if ( verbose ) cout << " Test 2) no fail here either " << endl ; 

    if (Errors < 0 ) {
      if (verbose ) cout << "Amesos_Klu" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestKlu relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestKlu relerror = " << relerror << endl ; 
      if (verbose) cout << " TestKlu maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestKlu maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestKlu NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
      exit( -13 ) ; 
    }
  }

#if 0
  //
  //     3)  ScaleMethod = 1 - argh I don't see how ScaleMEthod can work
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
    ParamList.set( "ScaleMethod", 1 );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 3) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Klu",
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
      if (verbose ) cout << "Amesos_Klu" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestKlu relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestKlu relerror = " << relerror << endl ; 
      if (verbose) cout << " TestKlu maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestKlu maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestKlu NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif
#if 0 

  This fails on Stratus  - see bug #1420

  //
  //     4)  ComputeTrueResidual==true
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
    ParamList.set( "ComputeTrueResidual", true );
    //    Teuchos::ParameterList& KluParams = ParamList.sublist("Klu") ;
    //    KluParams.set( "grid_mb", 3 );

    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 4) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
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
      if (verbose ) cout << "Amesos_Klu" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestKlu relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestKlu relerror = " << relerror << endl ; 
      if (verbose) cout << " TestKlu maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestKlu maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestKlu NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif

    return NumErrors; 
}

