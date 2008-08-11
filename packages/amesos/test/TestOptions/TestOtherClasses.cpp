#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "TestOtherClasses.h"
#include "PerformOneSolveAndTest.h"

int TestOtherClasses( const char* AmesosClass,
		      int EpetraMatrixType,
		      Epetra_CrsMatrix *& Amat, 
		      const bool transpose, 
		      const bool verbose, 
		      const int Levels,
		      const double Rcond,
		      bool RowMapEqualsColMap, 
		      bool TestAddZeroToDiag,
		      int ExpectedError,
		      double &maxrelerror, 
		      double &maxrelresidual,
		      int &NumTests ) {


  int iam = Amat->Comm().MyPID() ;  
  int NumErrors = 0 ;
  maxrelerror = 0.0;
  maxrelresidual = 0.0;
  const Epetra_Comm& Comm = Amat->Comm();

  bool MyVerbose = false ; // if set equal to verbose, we exceed the test harness 1 Megabyte limit
  string StringAmesosClass = AmesosClass ; 
  if ( AmesosClass ) MyVerbose = verbose ;    // Turn this on temporarily to debug Mumps on atlantis
  {
    Teuchos::ParameterList ParamList ;

    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", false );
    ParamList.set( "MaxProcs", 100000 );
    //  ParamList.print( cerr, 10 ) ; 

    double relerror;
    double relresidual;
   
    if (MyVerbose) cout << __FILE__ << "::" << __LINE__ << " AmesosClass= " << AmesosClass 
		      << " ParamList = " << ParamList 
		      << " transpose = " << transpose 
		      << " Levels = " << Levels 
		      << endl ; 

    int Errors = PerformOneSolveAndTest(AmesosClass,
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual, ExpectedError ) ; 

    if (MyVerbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
      if ( MyVerbose ) {
	cout << AmesosClass << " failed with error code " << Errors << " " << __FILE__ << "::" << __LINE__ << endl ; 
      }
    } else { 
      NumErrors += Errors ; 

      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( MyVerbose && Errors ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  string AC = AmesosClass ; 

  bool taucs = ( AC ==  "Amesos_Taucs" );
  bool klu = ( AC ==  "Amesos_Klu" );
  bool paraklete = ( AC ==  "Amesos_Paraklete" );
  bool mumps = ( AC ==  "Amesos_Mumps" );
  bool scalapack = ( AC ==  "Amesos_Scalapack" ) ;
  bool lapack = ( AC ==  "Amesos_Lapack" );


  //
  //  Testing AddZeroToDiag and AddToDiag 
  //  When AddZeroToDiag is true, the value of AddToDiag is added to every diagonal element whether or not 
  //    that element exists in the structure of the matrix.
  //  When AddZeroToDiag is false, the value of AddToDiag is added only to those diagonal elements 
  //    which are structually non-zero.
  //  Support for these two flags varies
  //
  //
  //  klu, superludist and parakalete support AddToDiag with or without AddZeroToDiag 
  //  scalapack and lapack, being dense codes, support AddToDiag, but only when AddZeroToDiag is set 
  //
  //  pardiso does not support AddToDiag - bug #1993 
  bool supports_AddToDiag_with_AddZeroToDiag = ( klu || paraklete || scalapack || lapack ) ; 
  bool supports_AddToDiag_with_when_AddZeroTo_Diag_is_false = ( klu  || paraklete  || mumps || taucs || lapack ) ; 


  if ( RowMapEqualsColMap && supports_AddToDiag_with_AddZeroToDiag && TestAddZeroToDiag ) {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", true );
    ParamList.set( "AddToDiag", 1.3e2 );

    //  ParamList.print( cerr, 10 ) ; 

    double relerror;
    double relresidual;
   
    if (MyVerbose) cout << __FILE__ << "::" << __LINE__ << " AmesosClass= " << AmesosClass 
		      << " ParamList = " << ParamList 
		      << " transpose = " << transpose 
		      << " Levels = " << Levels 
		      << endl ; 

    int Errors = PerformOneSolveAndTest(AmesosClass,
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual, ExpectedError ) ; 


    if (MyVerbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
      if ( MyVerbose ) {
	cout  << __FILE__ << "::" << __LINE__ 
	  << AmesosClass << " failed with error code " << Errors << " " << __FILE__ << "::" << __LINE__ << endl ; 
      }
    } else { 
      NumErrors += Errors ; 

      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( MyVerbose && Errors ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }



  }
  if ( RowMapEqualsColMap && supports_AddToDiag_with_when_AddZeroTo_Diag_is_false ) {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddToDiag", 1e2 );

    //  ParamList.print( cerr, 10 ) ; 

    double relerror;
    double relresidual;
   
    if (MyVerbose) cout << __FILE__ << "::" << __LINE__ << " AmesosClass= " << AmesosClass 
		      << " ParamList = " << ParamList 
		      << " transpose = " << transpose 
		      << " Levels = " << Levels 
		      << endl ; 

    int Errors = PerformOneSolveAndTest(AmesosClass,
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual, ExpectedError ) ; 


    if (MyVerbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
      if ( MyVerbose ) {
	cout  << __FILE__ << "::" << __LINE__ 
	  << AmesosClass << " failed with error code " << Errors << " " << __FILE__ << "::" << __LINE__ << endl ; 
      }
    } else { 
      NumErrors += Errors ; 

      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( MyVerbose && Errors ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }



  }

  //
  //     2)  Refactorize = true 
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "Refactorize", true );
      
    double relerror;
    double relresidual;
      
    if (MyVerbose) cout << __FILE__ << "::" << __LINE__ << " AmesosClass= " << AmesosClass 
		      << " ParamList = " << ParamList 
		      << " transpose = " << transpose 
		      << " Levels = " << Levels 
		      << endl ; 

    int Errors = PerformOneSolveAndTest(AmesosClass,
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual, ExpectedError ) ; 
      
    if (MyVerbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if (Errors < 0 ) {
      if (MyVerbose ) cout << AmesosClass << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( MyVerbose && Errors ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  //
  //     5)  MaxProcs = 2 
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "MaxProcs", 2 );   // Only affects Paraklete (maybe Mumps) also Superludist byt that isn't tested here 
      
    double relerror;
    double relresidual;
      
    if (MyVerbose) cout << __FILE__ << "::" << __LINE__ << " AmesosClass= " << AmesosClass 
		      << " ParamList = " << ParamList 
		      << " transpose = " << transpose 
		      << " Levels = " << Levels 
		      << endl ; 

    int Errors = PerformOneSolveAndTest(AmesosClass,
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual, ExpectedError ) ; 
      
    if (MyVerbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if (Errors < 0 ) {
      if (MyVerbose ) cout << AmesosClass << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( MyVerbose && Errors ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
  //
  //  ComputeTrueResidual is, by design, not quiet - it prints out the residual 
  //
#if 0
  //
  //     4)  ComputeTrueResidual==true
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "ComputeTrueResidual", true );
      
    double relerror;
    double relresidual;
      
    int Errors = PerformOneSolveAndTest(AmesosClass,
					EpetraMatrixType,
					Comm, 
					transpose, 
					MyVerbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual, ExpectedError ) ;

    if (Errors < 0 ) {
      if (MyVerbose ) cout << AmesosClass << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (MyVerbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( MyVerbose && Errors > 0 ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif


  return NumErrors; 
  }
