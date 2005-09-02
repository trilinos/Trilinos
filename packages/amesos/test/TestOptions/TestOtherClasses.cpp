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
		      double &maxrelerror, 
		      double &maxrelresidual,
		      int &NumTests ) {

  int iam = Amat->Comm().MyPID() ;  
  int NumErrors = 0 ;
  maxrelerror = 0.0;
  maxrelresidual = 0.0;
  const Epetra_Comm& Comm = Amat->Comm();

  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", true );
    ParamList.set( "MaxProcs", 100000 );
    //  ParamList.print( cerr, 10 ) ; 

    double relerror;
    double relresidual;
   
    if (verbose) cout << __FILE__ << "::" << __LINE__ << " AmesosClass= " << AmesosClass 
		      << " ParamList = " << ParamList 
		      << " transpose = " << transpose 
		      << " Levels = " << Levels 
		      << endl ; 

    int Errors = PerformOneSolveAndTest(AmesosClass,
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

    if (verbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
      if ( verbose ) {
	cout << AmesosClass << " failed with error code " << Errors << " " << __FILE__ << "::" << __LINE__ << endl ; 
      }
    } else { 
      NumErrors += Errors ; 

      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (verbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  if ( RowMapEqualsColMap ) {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddToDiag", 1e-3 );

    //  ParamList.print( cerr, 10 ) ; 

    double relerror;
    double relresidual;
   
    int Errors = PerformOneSolveAndTest(AmesosClass,
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


    if (verbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
      if ( verbose ) {
	cout << AmesosClass << " failed with error code " << Errors << " " << __FILE__ << "::" << __LINE__ << endl ; 
      }
    } else { 
      NumErrors += Errors ; 

      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (verbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors ) {
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
      
    int Errors = PerformOneSolveAndTest(AmesosClass,
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
      
    if (verbose  || ( Errors && iam==0 )  ) cout << __FILE__ << "::" << __LINE__ 
				  << " AmesosClass= " << AmesosClass 
				  << " Errors = " << Errors 
				  << endl ; 

    if (Errors < 0 ) {
      if (verbose ) cout << AmesosClass << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

    }
    if (verbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors ) {
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

    }
    if (verbose)  cout << " TestOtherClasses " << AmesosClass << "" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif


  return NumErrors; 
  }
