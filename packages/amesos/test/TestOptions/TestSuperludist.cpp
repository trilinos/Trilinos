#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "PerformOneSolveAndTest.h"


//
//  Returns the number of failures.
//  Note:  If AMESOS_SUPERLUDIST is not supported, TestSuperludist() will 
//  always return 0
//
//  TestSuperludist performs the following tests:
//                         Redistribute   AddZeroToDiag   SUB:  ReuseSymbolic MaxProcesses
//                            true           true                   true           2
//                            true           true                   false          2
//                            true           true                   false          2
//                            true           true                   true           1
//                            true           true                   false          1
//                            true           true                   false          2
//                            true           false                  true           1
//                            true           false                  true           2
//                            true           false                  false          1
//                            true           false                  false          2
//                            false/true     true                   true           1
//                            false          true                   true           2
//                            false          true                   false          1
//                            false          true                   false          2
//                            false          false                  true           1
//                            false          false                  true           2
//                            false          false                  false          1
//                            false          false                  false          2
//   


int TestSuperludist( Epetra_CrsMatrix *& Amat, 
		     int EpetraMatrixType,
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
  double relerror;
  double relresidual;
   

  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
    
    int Errors = PerformOneSolveAndTest("Amesos_Superludist",
					EpetraMatrixType,
					Comm, 
					transpose, 
					verbose,
					ParamList, 
					Amat, 
					Levels,
					Rcond, 
					relerror, 
					relresidual) ;

    if ( Errors < 0 ) {
      NumErrors++;
      NumTests++ ; 
      if ( verbose ) {
	cout << "Amesos_Superludist failed with error code " << Errors<< endl ; 
      }
    } else { 

      NumErrors += Errors ; 
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
					  
					  
    }
      
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
	
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
    }
      
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  

    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 1 );
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
    {
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
  
    {
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
      maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
      maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
      NumTests++ ; 
	
      //      NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
    }
  
    return NumErrors; 
  }

}
