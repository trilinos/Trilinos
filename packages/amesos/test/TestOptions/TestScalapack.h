#include "Epetra_CrsMatrix.h"
 
int TestScalapack( Epetra_CrsMatrix *& Amat, 
		     int EpetraMatrixType,
		     bool transpose, 
		     bool verbose, 
		     int Levels,
		     const double Rcond,
		     double &maxrelerror, 
		     double &maxrelresidual,
		     int &NumTests) ;


#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif
