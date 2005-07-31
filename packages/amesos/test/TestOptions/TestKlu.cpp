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
//     2A)  AddToDiag = 1e-12
#if 0
//     3)  ScaleMethod = 1 - argh I don't see how ScaleMEthod can work
//     4)  ComputeTrueResidual==true
#endif
//
int TestKlu( Epetra_CrsMatrix *& Amat, 
	     int EpetraMatrixType,
	     const bool transpose, 
	     const bool verbose, 
	     const int Levels,
	     const double Rcond,
	     Teuchos::ParameterList ParamList,
	     bool RowMapEqualsColMap, 
	     double &maxrelerror, 
	     double &maxrelresidual,
	     int &NumTests ) {
  
  int NumErrors = 0 ;
  maxrelerror = 0.0;
  maxrelresidual = 0.0;
  const Epetra_Comm& Comm = Amat->Comm();

  if ( verbose ) ParamList.set( "DebugLevel", 1 );
  if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
  else ParamList.set( "OutputLevel", 1 );

  //
  //     1)  no parameters

  {
    Teuchos::ParameterList InternalParamList = ParamList ; 
      
    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual ) ;

      
    if ( Amat->Comm().MyPID() == 0 && Errors ) {
      cout << __FILE__ << "::"  << __LINE__ 
	   << "Amesos_Klu failed with error code " << Errors<< endl ; 
      }
    if (Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
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
    if ( Amat->Comm().MyPID() == 0 && Errors > 0 ) {
      cout << "Amesos_Klu" 
	   << __FILE__ << "::"  << __LINE__ 
	   << " Errors = " <<  Errors 
	   << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  //
  //     2)  Refactorize = true 
  {
    Teuchos::ParameterList InternalParamList = ParamList ; 

    InternalParamList.set( "Refactorize", true );
      
    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual ) ; 
      
    if (  Amat->Comm().MyPID() == 0 && Errors ) {
      cout << __FILE__ << "::"  << __LINE__ 
	   << "Amesos_Klu failed with error code " << Errors<< endl ; 
      }
    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
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
    if (  Amat->Comm().MyPID() == 0 && Errors > 0 ) {
      cout << "Amesos_Klu" 
	   << __FILE__ << "::"  << __LINE__ 
	   << " Errors = " <<  Errors 
	   << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
  //
  //     2a) TrustMe = true 
  //         Note:  Requires Optimized Storage (i.e. EpetraMatrixType == 2 ) 
  //                and does not support reindexing
  bool ReIndex = ParamList.get( "Reindex", false );
  bool DontTrustMe = ParamList.get( "DontTrustMe", false );
  if ( EpetraMatrixType == 2 && ! ReIndex && ! DontTrustMe ) {
    Teuchos::ParameterList InternalParamList = ParamList ; 

    InternalParamList.set( "TrustMe", true );
      
    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual ) ; 
      
    if (  Amat->Comm().MyPID() == 0 && Errors ) {
      cout << __FILE__ << "::"  << __LINE__ 
	   << "Amesos_Klu failed with error code " << Errors<< endl ; 
      }
    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
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
    if (  Amat->Comm().MyPID() == 0 && Errors > 0 ) {
      cout << "Amesos_Klu" 
	   << __FILE__ << "::"  << __LINE__ 
	   << " Errors = " <<  Errors 
	   << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
  //
  //     2A)  AddToDiag 
  if (RowMapEqualsColMap ) {
    Teuchos::ParameterList InternalParamList = ParamList ; 
    InternalParamList.set( "Refactorize", true );
    InternalParamList.set( "AddToDiag", 1e-2 );

    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual ) ; 
      
    if (  Amat->Comm().MyPID() == 0 && Errors ) {
      cout << __FILE__ << "::"  << __LINE__ 
	   << "Amesos_Klu failed with error code " << Errors<< endl ; 
      }
    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
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
    if ( Comm.MyPID() == 0 && Errors > 0 ) {
      cout << "Amesos_Klu" 
	   << __FILE__ << "::"  << __LINE__ 
	   << " Errors = " <<  Errors 
	   << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

#if 0
  //
  //     3)  ScaleMethod = 1 - argh I don't see how ScaleMEthod can work
  {
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
    ParamList.set( "ScaleMethod", 1 );
      
    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					Comm, 
					transpose, 
					verbose,
					InternalParamList, 
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
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    else ParamList.set( "OutputLevel", 1 );
    ParamList.set( "ComputeTrueResidual", true );
    //    Teuchos::ParameterList& KluParams = ParamList.sublist("Klu") ;
    //    KluParams.set( "grid_mb", 3 );

    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					InternalParamList, 
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

