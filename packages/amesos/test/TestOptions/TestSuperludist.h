#include "Epetra_CrsMatrix.h"

int TestSuperludist( Epetra_CrsMatrix *& Amat,
                     int EpetraMatrixType,
                     bool transpose,
                     bool verbose,
                     int Levels,
                     const double Rcond,
                     double &maxrelerror,
                     double &maxrelresidual,
                     const char *filename,
                     int &NumTests) ;

