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
//  Note:  If AMESOS_SUPERLUDIST is not supported, PerformOneSolveAndTest() will 
//  always return 0
//
//  Still have to decide where we are going to check the residual.  
//
//  The following table shows the variable names that we use for 
//  each of the three phases:  
//     compute - which computes the correct value of b
//     solve - which solves for x in  A' A' A x = b 
//     check - which computes bcheck = A' A' A x 
//
//  For ill-conditioned matrices we restrict the test to one or two 
//  solves, by setting Levels to 1 or 2 on input to this routine.
//  When Levels is less than 3, some of the transformations
//  shown in the table as "->" and "<-" are not performed, instead 
//  a direct copy is made.
//
//  In the absence of roundoff, each item in a given column should 
//  be identical.  
//
//  If Levels = 3, we compute and solve A' A' A x = b and hence all 
//  of the transformations are performed
//
//  If Levels = 2, the transformations shown in the first column of 
//  transformations (labelled Levels>=3) are replaced by a direct copy.
//
//  If Levels = 1, only the transformations shown in the third column
//  are performed, the others being replaced by direct copies.
//  
//                           Levels>=3    Levels>=2
//                              A'         A'            A
//  compute             xexact  ->  cAx    ->     cAAx   ->       b 
//  solve               x       <-  sAx    <-     sAAx   <-       b
//  check               x       ->  kAx    ->     kAAx   ->  bcheck
//
//  Note that since Levels 2 and 3 use the same A, there is no need to 
//  call NumericFactorization() between the second and third call to Solve. 
//   


int TestSuperludist( Epetra_CrsMatrix *& Amat, 
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
    
    NumErrors += PerformOneSolveAndTest("Amesos_Superludist",
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
    
    if (verbose) cout << " TestSuperludist relresidual = " <<relresidual << endl ; 
    if (verbose) cout << " TestSuperludist relerror = " << relerror << endl ; 
    if (verbose) cout << " TestSuperludist maxrelresidual = " << maxrelresidual << endl ; 
    if (verbose) cout << " TestSuperludist maxrelerror = " << maxrelerror << endl ; 
    
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

#if 0
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", true );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  

  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", true );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", true );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 1 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
  
  
  {
    Teuchos::ParameterList ParamList ;
    ParamList.set( "Redistribute", false );
    ParamList.set( "AddZeroToDiag", false );
    Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
    SuperludistParams.set( "ReuseSymbolic", false );
    SuperludistParams.set( "MaxProcesses", 2 );
    //  ParamList.print( cerr, 10 ) ; 
   
    NumErrors += PerformOneSolveAndTest( Comm, ParamList ) ; 
  }
#endif
  
  return NumErrors; 
}
