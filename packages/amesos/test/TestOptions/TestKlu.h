#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"
 
int TestKlu( Epetra_CrsMatrix *& Amat, 
	     int EpetraMatrixType,
	     const bool transpose, 
	     const bool verbose, 
	     const int Levels,
	     const double Rcond,
	     Teuchos::ParameterList ParamList, 
	     bool RowMapEqualsColMap, 
             bool TestAddZeroToDiag,
	     int ExpectedError,
	     double &maxrelerror, 
	     double &maxrelresidual,
	     int &NumTests ) ;

