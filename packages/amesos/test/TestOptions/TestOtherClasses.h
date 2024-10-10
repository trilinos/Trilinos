#include "Epetra_CrsMatrix.h"
 
int TestOtherClasses( const char* AmesosClass,
		     int EpetraMatrixType,
		      Epetra_CrsMatrix *& Amat, 
		      const bool transpose, 
		      const bool verbose, 
		      const int Levels,
		      const double Rcond,
		      bool RowMapEqualsColMap, 
		      bool TestAddZeroToDiag,
		      int ExpectedError,
		      double &maxrelerror, 
		      double &maxrelresidual,
		      int &NumTests ) ;


#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif
