//
//  Knonw Bugs: 
//     Fails in test #6 
//     Tested only on SuperLU.rua
//     Monstrously too many prints
//     Fails on an assert when teseint glarger matrix
//        Also in test #6 - try without test #6 
//
//     Epetra_CHK_ERRs look like negative values on return - and hence
//     look like Amesos_Scalapack is not built in. 
//        PROPOSED SOLUTION:  Set a BUILT  flag for each AmesosClass 
//     verbose is not controlling print outs from Amesos_Scalapack -
//       they are always on  -DONE
//     Failue in PDGETRS - FIXED 
//     debug is not being turned off in Amesos_ScaLAPACK - FIXED
//
//
//








#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "PerformOneSolveAndTest.h"
#include "TestScalapack.h"

//
//  Tests:
//     1)  MaxProcs==100000
//     2)  ComputeVectorNorms==true,  MaxProcs==2
//     3)  ComputeVectorNorms==true,  MaxProcs==2, grid_mb_ =2, grid_nb_=2
//     4)  ComputeTrueResidual==true,  MaxProcs==1000 grid_mb_ =3, grid_nb_=3
//     5)  "2D distribution"=false, ComputeTrueResidual=true, ComputeVectorNorms=true
//     6)  no parameters
//
int TestScalapack( Epetra_CrsMatrix *& Amat, 
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

#if 1
  //
  //     1)  MaxProcs==100000

  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    ParamList.set( "MaxProcs", 100000 );
    //  ParamList.print( cerr, 10 ) ; 
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 1) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Scalapack",
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
	cout << "Amesos_Scalapack failed with error code " << Errors<< endl ; 
      }
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestScalapack relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestScalapack relerror = " << relerror << endl ; 
      if (verbose) cout << " TestScalapack maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestScalapack maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestScalapack NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Scalapack" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  //
  //     2)  ComputeVectorNorms==true,  MaxProcs==2
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    ParamList.set( "MaxProcs", 2 );
    ParamList.set( "ComputeTrueResidual", true );
    //  ParamList.print( cerr, 10 ) ; 
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 2) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Scalapack",
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
      if (verbose ) cout << "Amesos_Scalapack" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestScalapack relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestScalapack relerror = " << relerror << endl ; 
      if (verbose) cout << " TestScalapack maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestScalapack maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestScalapack NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Scalapack" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  //
  //     3)  ComputeVectorNorms==true,  MaxProcs==2, grid_mb_ =2, grid_nb_=2
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    ParamList.set( "MaxProcs", 2 );
    ParamList.set( "ComputeTrueResidual", true );
    //  ParamList.print( cerr, 10 ) ; 
    Teuchos::ParameterList& ScalapackParams = ParamList.sublist("Scalapack") ;
    ScalapackParams.set( "grid_mb", 2 );
    ScalapackParams.set( "grid_nb", 2 );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 3) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Scalapack",
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
      if (verbose ) cout << "Amesos_Scalapack" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestScalapack relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestScalapack relerror = " << relerror << endl ; 
      if (verbose) cout << " TestScalapack maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestScalapack maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestScalapack NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Scalapack" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  //
  //     4)  ComputeTrueResidual==true,  MaxProcs==1000 grid_mb_ =3, grid_nb_=3
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    ParamList.set( "MaxProcs", 1000 );
    ParamList.set( "ComputeTrueResidual", true );
    Teuchos::ParameterList& ScalapackParams = ParamList.sublist("Scalapack") ;
    ScalapackParams.set( "grid_mb", 3 );
    ScalapackParams.set( "grid_nb", 3 );

    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 4) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Scalapack",
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
      if (verbose ) cout << "Amesos_Scalapack" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestScalapack relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestScalapack relerror = " << relerror << endl ; 
      if (verbose) cout << " TestScalapack maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestScalapack maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestScalapack NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Scalapack" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

  //
  //     5)  "2D distribution"=false, ComputeTrueResidual=true, ComputeVectorNorms=true
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
    ParamList.set( "MaxProcs", 1000 );
    ParamList.set( "ComputeTrueResidual", true );
    Teuchos::ParameterList& ScalapackParams = ParamList.sublist("Scalapack") ;
    ScalapackParams.set( "2D distribution", false );

      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 5) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Scalapack",
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
      if (verbose ) cout << "Amesos_Scalapack" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestScalapack relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestScalapack relerror = " << relerror << endl ; 
      if (verbose) cout << " TestScalapack maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestScalapack maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestScalapack NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Scalapack" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }

#endif
#if 1
  //
  //     6)  no parameters
  {
    Teuchos::ParameterList ParamList ;
    if ( verbose ) ParamList.set( "DebugLevel", 1 );
    if ( ! verbose ) ParamList.set( "OutputLevel", 0 );
      
    double relerror;
    double relresidual;
      
    if ( verbose ) cout << " Test 6) no fail yet " << endl ; 

    int Errors = PerformOneSolveAndTest("Amesos_Scalapack",
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
      if (verbose ) cout << "Amesos_Scalapack" << " not built in this executable " << endl ; 
      return 0 ; 
    } else { 
      NumErrors += Errors ; 
	
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 

      if (verbose) cout << " TestScalapack relresidual = " <<relresidual << endl ; 
      if (verbose) cout << " TestScalapack relerror = " << relerror << endl ; 
      if (verbose) cout << " TestScalapack maxrelresidual = " << maxrelresidual << endl ; 
      if (verbose) cout << " TestScalapack maxrelerror = " << maxrelerror << endl ; 
	
    }
    if (verbose)  cout << " TestScalapack NumErrors = " << NumErrors << endl ; 
    if ( verbose && Errors > 0 ) {
      cout << "Amesos_Scalapack" << " failed with transpose = " << 
	(transpose?"true":"false") << endl ;  
    }
  }
#endif


    return NumErrors; 
}

