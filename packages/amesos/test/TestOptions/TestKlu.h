#include "Epetra_CrsMatrix.h"
 
int TestKlu( Epetra_CrsMatrix *& Amat, 
	     int EpetraMatrixType,
	     const bool transpose, 
	     const bool verbose, 
	     const int Levels,
	     const double Rcond,
	     bool RowMapEqualsColMap, 
	     double &maxrelerror, 
	     double &maxrelresidual,
	     int &NumTests ) ;

