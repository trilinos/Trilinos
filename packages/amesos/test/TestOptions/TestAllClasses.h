#include "Epetra_CrsMatrix.h"
#include <vector>

int TestAllClasses(vector<string> AmesosClasses,
		   vector<bool> AmesosClassesInstalled,
		   Epetra_CrsMatrix *& Amat, 
		   bool transpose, 
		   bool verbose, 
		   bool symmetric, 
		   int Levels,
		   const double Rcond,
		   double &maxrelerror, 
		   double &maxrelresidual,
		   int &NumTests) ;

