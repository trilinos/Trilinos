#include "Epetra_Comm.h"
#include "Amesos_Parameter_List.h"
#include "Amesos_Factory.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "TestOtherClasses.h"
#include "PerformOneSolveAndTest.h"

int TestOtherClasses(AmesosClassType AmesosClass,
		     Epetra_CrsMatrix *& Amat, 
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

{
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 

   double relerror;
   double relresidual;
   
   NumErrors += PerformOneSolveAndTest(AmesosClass,
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

   if (verbose) cout << " TestOtherClasses relresidual = " <<relresidual << endl ; 
   if (verbose) cout << " TestOtherClasses relerror = " << relerror << endl ; 
   if (verbose) cout << " TestOtherClasses maxrelresidual = " << maxrelresidual << endl ; 
   if (verbose) cout << " TestOtherClasses maxrelerror = " << maxrelerror << endl ; 

 }
  cout << " TestOtherClasses NumErrors = " << NumErrors << endl ; 
  
 return NumErrors; 
}
