#include "Epetra_CrsMatrix.h"
 
int TestSuperludist( Epetra_CrsMatrix *& Amat, 
		     bool transpose, 
		     bool verbose, 
		     int Levels,
		     const double Rcond,
		     double &maxrelerror, 
		     double &maxrelresidual,
		     int &NumTests) ;

