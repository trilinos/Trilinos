#include "Epetra_CrsMatrix.h"
 
int TestAllClasses(Epetra_CrsMatrix *& Amat, 
		   bool transpose, 
		   bool verbose, 
		   bool symmetric, 
		   int Levels,
		   const double Rcond,
		   double &maxrelerror, 
		   double &maxrelresidual,
		   int &NumTests) ;

