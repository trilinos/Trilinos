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
    if ( verbose ) ParamList.set( "DebugLevel", 1 );

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

      if (verbose) cout << " TestOtherClasses relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestOtherClasses relerror = " << relerror << endl ; 
      if (verbose) cout << " TestOtherClasses maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestOtherClasses maxrelerror = " << maxrelerror << endl ; 

    }
    if (verbose)  cout << " TestOtherClasses" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  if ( RowMapEqualsColMap ) {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddToDiag", 1e-3 );
    if ( verbose ) ParamList.set( "DebugLevel", 1 );

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

      if (verbose) cout << " TestOtherClasses relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestOtherClasses relerror = " << relerror << endl ; 
      if (verbose) cout << " TestOtherClasses maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestOtherClasses maxrelerror = " << maxrelerror << endl ; 

    }
    if (verbose)  cout << " TestOtherClasses" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }



  }
  //
  //     2)  Refactorize = true 
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    ParamList.set( "Refactorize", true );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 2) no fail yet " << endl ; 

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
      
    if ( verbose ) cout << " Test 2) no fail down here " << endl ; 

    if (Errors < 0 ) {
      if (verbose ) cout << AmesosClass << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << AmesosClass << "  relresidual = " <<relresidual << endl ; 
      if (verbose) cout << AmesosClass << "  relerror = " << relerror << endl ; 
      if (verbose) cout << AmesosClass << "  maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << AmesosClass << "  maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestOtherClasses" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
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
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    ParamList.set( "ComputeTrueResidual", true );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 2) no fail yet " << endl ; 

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

      if (verbose) cout << AmesosClass << "  relresidual = " <<relresidual << endl ; 
      if (verbose) cout << AmesosClass << "  relerror = " << relerror << endl ; 
      if (verbose) cout << AmesosClass << "  maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << AmesosClass << "  maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestOtherClasses" << "::" << __LINE__ << " NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << AmesosClass << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif


  return NumErrors; 
  }
