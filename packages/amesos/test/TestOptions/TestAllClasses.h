#include "Epetra_CrsMatrix.h"
#include <vector>

int TestAllClasses( const vector<string> AmesosClasses,
		    const vector<bool> AmesosClassesInstalled,
		    Epetra_CrsMatrix *& Amat, 
		    const bool transpose, 
		    const bool verbose, 
		    const bool symmetric, 
		    const int Levels,
		    const double Rcond,
		    double &maxrelerror, 
		    double &maxrelresidual,
		    int &NumTests) ;

