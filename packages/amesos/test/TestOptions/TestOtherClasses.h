#include "Epetra_CrsMatrix.h"
 
int TestOtherClasses( const char* AmesosClass,
		      Epetra_CrsMatrix *& Amat, 
		      const bool transpose, 
		      const bool verbose, 
		      const int Levels,
		      const double Rcond,
		      double &maxrelerror, 
		      double &maxrelresidual,
		      int &NumTests ) ;

