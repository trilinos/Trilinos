#include "Epetra_Comm.h"
#include "Amesos_Parameter_List.h"
#include "Amesos_Factory.h"
#include "Epetra_CrsMatrix.h"

int PerformOneSolveAndTest(AmesosClassType AmesosClass,
			   const Epetra_Comm &Comm, 
			   bool transpose, 
			   bool verbose, 
			   AMESOS::Parameter::List ParamList, 
			   Epetra_CrsMatrix *& Amat, 
			   int Levels, 
			   const double Rcond,
			   double& relerror,
			   double& relresidual) ;
