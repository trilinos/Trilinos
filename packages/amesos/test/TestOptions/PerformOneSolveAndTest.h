#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos_Factory.h"
#include "Epetra_CrsMatrix.h"

int PerformOneSolveAndTest(AmesosClassType AmesosClass,
			   const Epetra_Comm &Comm, 
			   bool transpose, 
			   bool verbose, 
			   Teuchos::ParameterList ParamList, 
			   Epetra_CrsMatrix *& Amat, 
			   int Levels, 
			   const double Rcond,
			   double& relerror,
			   double& relresidual) ;
