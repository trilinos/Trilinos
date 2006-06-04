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
             bool TestAddZeroToDiag,
	     int ExpectedError,
	     double &maxrelerror, 
	     double &maxrelresidual,
	     int &NumTests ) {
  
  bool MyVerbose = false ; //  if set equal to verbose, we exceed thee test harness 1 Megabyte limit

  int NumErrors = 0 ;
  maxrelerror = 0.0;
  maxrelresidual = 0.0;
  const Epetra_Comm& Comm = Amat->Comm();

  //
  //     1)  no parameters

  {
    Teuchos::ParameterList InternalParamList = ParamList ; 
      
    double relerror;
    double relresidual;
      
    if ( MyVerbose ) cout  << __FILE__ << "::"  << __LINE__ 
      << " InternalParamList = " <<
		     InternalParamList <<  endl ; 
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual,
					ExpectedError ) ;

      
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

    }
    if (MyVerbose)  cout << " TestKlu NumErrors = " 
		       << NumErrors << " "
		       << __FILE__ << "::" << __LINE__ 
		       << endl ; 
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
    if ( MyVerbose ) cout  << __FILE__ << "::"  << __LINE__ 
			 << " InternalParamList = " <<
		     InternalParamList <<  endl ; 
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual,
					ExpectedError ) ;
      
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

    }
    if (MyVerbose)  cout << " TestKlu NumErrors = " 
		       << NumErrors << " "
		       << __FILE__ << "::" << __LINE__ 
		       << endl ; 
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

    if ( MyVerbose ) cout  << __FILE__ << "::"  << __LINE__ 
      << " InternalParamList = " <<
		     InternalParamList <<  endl ; 
      
    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual,
					ExpectedError ) ;
 
      
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

    }
    if (MyVerbose)  cout << " TestKlu NumErrors = " 
		       << NumErrors << " "
		       << __FILE__ << "::" << __LINE__ 
		       << endl ; 
    if (  Amat->Comm().MyPID() == 0 && Errors > 0 ) {
      cout << "Amesos_Klu" 
	   << __FILE__ << "::"  << __LINE__ 
	   << " Errors = " <<  Errors 
	   << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
  

  if ( ExpectedError == 0 ) {
    //
    //     2A)  AddToDiag 
    if (RowMapEqualsColMap ) {
      Teuchos::ParameterList InternalParamList = ParamList ; 
      InternalParamList.set( "Refactorize", true );
      //    InternalParamList.set( "AddZeroToDiag", true );
      InternalParamList.set( "AddToDiag", 1e2 );
      
      double relerror;
      double relresidual;
      
      if ( MyVerbose ) cout  << __FILE__ << "::"  << __LINE__ 
			     << " InternalParamList = " <<
			 InternalParamList <<  endl ; 
      
      int Errors = PerformOneSolveAndTest("Amesos_Klu",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
					  InternalParamList, 
					  Amat, 
					  Levels,
					  Rcond, 
					  relerror, 
					  relresidual,
					  ExpectedError ) ;
      
      
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
	
      }
      if (MyVerbose)  cout << " TestKlu NumErrors = " 
			   << NumErrors << " "
			   << __FILE__ << "::" << __LINE__ 
			   << endl ; 
      if ( Comm.MyPID() == 0 && Errors > 0 ) {
	cout << "Amesos_Klu" 
	     << __FILE__ << "::"  << __LINE__ 
	     << " Errors = " <<  Errors 
	     << " failed with transpose = " << 
	  (transpose?"true":"false") << endl ;  
      }
    }

    //
    //     2B)  AddToDiag with AddZeroToDiag
    if (RowMapEqualsColMap && TestAddZeroToDiag ) {
      Teuchos::ParameterList InternalParamList = ParamList ; 
      InternalParamList.set( "Refactorize", true );
      InternalParamList.set( "AddZeroToDiag", true );
      InternalParamList.set( "AddToDiag", 1e2 );
      
      double relerror;
      double relresidual;
      
      if ( MyVerbose ) cout  << __FILE__ << "::"  << __LINE__ 
			     << " InternalParamList = " <<
			 InternalParamList <<  endl ; 
      
      int Errors = PerformOneSolveAndTest("Amesos_Klu",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
					  InternalParamList, 
					  Amat, 
					  Levels,
					  Rcond, 
					  relerror, 
					  relresidual,
					  ExpectedError ) ;
      
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
	
      }
      if (MyVerbose)  cout << " TestKlu NumErrors = " 
			   << NumErrors << " "
			   << __FILE__ << "::" << __LINE__ 
			   << endl ; 
      if ( Comm.MyPID() == 0 && Errors > 0 ) {
	cout << "Amesos_Klu" 
	     << __FILE__ << "::"  << __LINE__ 
	     << " Errors = " <<  Errors 
	     << " failed with transpose = " << 
	  (transpose?"true":"false") << endl ;  
      }
    }
    
    //
    //     2C)   AddZeroToDiag without AddToDiag 
    if (RowMapEqualsColMap  ) {
      Teuchos::ParameterList InternalParamList = ParamList ; 
      InternalParamList.set( "Refactorize", true );
      InternalParamList.set( "AddZeroToDiag", true );

      double relerror;
      double relresidual;
      
      if ( MyVerbose ) cout  << __FILE__ << "::"  << __LINE__ 
			     << " InternalParamList = " <<
			 InternalParamList <<  endl ; 
      
      int Errors = PerformOneSolveAndTest("Amesos_Klu",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
					  InternalParamList, 
					  Amat, 
					  Levels,
					  Rcond, 
					  relerror, 
					  relresidual,
					  ExpectedError ) ;
      
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
	
      }
      if (MyVerbose)  cout << " TestKlu NumErrors = " 
			   << NumErrors << " "
			   << __FILE__ << "::" << __LINE__ 
			   << endl ; 
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
    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					Comm, 
					transpose, 
					MyVerbose,
					InternalParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual,
					ExpectedError ) ;
 
      
    if (Errors < 0 ) {
      if (MyVerbose ) cout << "Amesos_Klu" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestKlu NumErrors = " 
		       << NumErrors << " "
		       << __FILE__ << "::" << __LINE__ 
		       << endl ; 
    if ( MyVerbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif
#if 0

  //  This fails on Stratus  - see bug #1420 - now marked as a duplicate of bug 1417

  //
  //     4)  ComputeTrueResidual==true
  {
    ParamList.set( "ComputeTrueResidual", true );
    //    Teuchos::ParameterList& KluParams = ParamList.sublist("Klu") ;
    //    KluParams.set( "grid_mb", 3 );

    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest("Amesos_Klu",
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual,
					ExpectedError ) ;

      
    if (Errors < 0 ) {
      if (MyVerbose ) cout << "Amesos_Klu" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestKlu NumErrors = " << NumErrors 
		       << " " << __FILE__ << "::" << __LINE__ << endl ; 
    if ( MyVerbose && Errors > 0 ) {
      cout << "Amesos_Klu" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif
  }
    return NumErrors; 
}

