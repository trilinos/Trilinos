#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"

int PerformOneSolveAndTest(const char* AmesosClass,
			   int EpetraMatrixType,
			   const Epetra_Comm &Comm, 
			   bool transpose, 
			   bool verbose, 
			   Teuchos::ParameterList ParamList, 
			   Epetra_CrsMatrix *& Amat, 
			   int Levels, 
			   const double Rcond,
			   double& relerror,
			   double& relresidual,
			   int ExpectedError = 0) ;

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif
