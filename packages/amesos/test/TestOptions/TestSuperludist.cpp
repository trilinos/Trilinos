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
//  Test number:
//   1 disabled AddToDiag=100 true           true                   true           2
//   2 disabled AddToDiag=100 true           false                  true           2
//   3                        true           true                   false          2
//   4                        true           true                   true           1
//   5                        true           true                   false          1
//   6                        true           true                   false          10
//   7                        true           false                  true           -1
//   8                        true           false                  true           -2
//   9                        true           false                  true           -3
//   10                       true           false                  false          4  
//   11                       false/true     true                   true           4
//   12                       false/true     true                   true           2
//     Test #12 appears to duplicate test #11 and perhaps #4 
//   13                       false/true     true                   false          1
//   14                       false/true     true                   false          2
//   15                       false/true     false                  true           1
//   16                       false/true     false                  true           2
//   17 SamePattern           true           false                  true           10
//   18 RowPerm - NATURAL     true           false                  false          10
//   19 RowPerm - LargeDiag   true           false                  false          10
//   20 RowPerm - NATURAL     true           false                  false          10
//   21 RowPerm - LargeDiag   true           false                  false          10
//   22 RowPerm - TinyPivot=t true           false                  false          10
//   23 RowPerm - TinyPivot=f true           false                  false          10
//   


int TestSuperludist( Epetra_CrsMatrix *& Amat, 
		     int EpetraMatrixType,
		     bool transpose, 
		     bool verbose, 
		     int Levels,
		     const double Rcond,
		     double &maxrelerror, 
		     double &maxrelresidual,
		     char *filename,
		     int &NumTests ) {

  std::string StringFilename = filename ; 
  bool ImpcolB = ( StringFilename.find("ImpcolB") < StringFilename.find("xdz_notaname_garbage") );
  int NumErrors = 0 ;
  maxrelerror = 0.0;
  maxrelresidual = 0.0;
  const Epetra_Comm& Comm = Amat->Comm();
  double relerror;
  double relresidual;
   

  {
    bool MyVerbose = false ; // if set to verbose - we exceed the test harness 1 Megabyte limit
    //  bool MyVerbose = verbose ; // if set to verbose - we exceed the test harness 1 Megabyte limit

    //
    //  Bug #1990 - AddToDiag fails in Amesos_Superludist 
    //
#if 0 
    //  Test #1 - disabled - bug #1990
    {
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      ParamList.set( "AddToDiag", 1e2 );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( std::cerr, 10 ) ; 
      
      const int ExpectedError = 0 ;
      int Errors = PerformOneSolveAndTest("Amesos_Superludist",
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
      
      if ( Errors < 0 ) {
	NumErrors++;
	NumTests++ ; 
	if ( MyVerbose ) {
	  std::cout << "Amesos_Superludist failed with error code " << Errors<< std::endl ; 
	}
      } else { 
	
	NumErrors += Errors ; 
	maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
	maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
	NumTests++ ; 
	
					  
      }
    }
    {

      // Test #2 - disabled -  bug #1990 - AddToDiag fails 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      ParamList.set( "AddToDiag", 1e2 );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( std::cerr, 10 ) ; 
      
      int Errors = PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
					  ParamList, 
					  Amat, 
					  Levels,
					  Rcond, 
					  relerror, 
					  relresidual) ;
      
      if ( Errors < 0 ) {
	NumErrors++;
	NumTests++ ; 
	if ( MyVerbose ) {
	  std::cout << "Amesos_Superludist failed with error code " << Errors<< std::endl ; 
	}
      } else { 
	
	NumErrors += Errors ; 
	maxrelerror = EPETRA_MAX( relerror, maxrelerror ) ; 
	maxrelresidual = EPETRA_MAX( relresidual, maxrelresidual ) ; 
	NumTests++ ; 
	
					  
      }
    }
#endif     
    {
      // test #3 - 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( std::cerr, 10 ) ; 
	
    if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
      << " ParamList = " <<
		     ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  test #4 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  test #5
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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

      //  test #6 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 10 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #7 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", -1 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #8 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", -2 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #9 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", -3 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #10
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 4 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #11 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 4 );
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #12 
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 4 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #13 
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      // Test #14 
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", true );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #15 
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 1 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #16 
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #17 
      Teuchos::ParameterList ParamList ;
      if ( Amat->RowMatrixRowMap().LinearMap() == false )   // bug #1408
	ParamList.set( "Redistribute", true );
      else
	ParamList.set( "Redistribute", false );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "Fact", "SamePattern" );
      SuperludistParams.set( "MaxProcesses", 2 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
  
    if (!ImpcolB )   // ImpcolB fails if the NATURAL order - i.e. no pivoting - is chosen  
    {
      //  Test #18 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "RowPerm", "NATURAL" );
      SuperludistParams.set( "MaxProcesses", 10 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #19 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", false );
      SuperludistParams.set( "RowPerm", "LargeDiag" );
      SuperludistParams.set( "MaxProcesses", 10 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
  
  
    if (!ImpcolB )   // ImpcolB fails if the NATURAL order - i.e. no pivoting - is chosen  
    {
      //  Test #20 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "RowPerm", "NATURAL" );
      SuperludistParams.set( "MaxProcesses", 10 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #21
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "RowPerm", "LargeDiag" );
      SuperludistParams.set( "MaxProcesses", 10 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #22
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "ReplaceTinyPivot", true );
      SuperludistParams.set( "MaxProcesses", 10 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
      //  Test #23
      Teuchos::ParameterList ParamList ;
      ParamList.set( "Redistribute", true );
      ParamList.set( "AddZeroToDiag", false );
      Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
      SuperludistParams.set( "ReuseSymbolic", true );
      SuperludistParams.set( "ReplaceTinyPivot", false );
      SuperludistParams.set( "MaxProcesses", 10 );
      //  ParamList.print( std::cerr, 10 ) ; 
   
      if ( MyVerbose ) std::cout  << __FILE__ << "::"  << __LINE__ 
			   << " ParamList = " <<
		       ParamList <<  std::endl ; 
      
      NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
					  EpetraMatrixType,
					  Comm, 
					  transpose, 
					  MyVerbose,
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
