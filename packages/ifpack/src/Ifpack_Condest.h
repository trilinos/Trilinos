#ifndef IFPACK_CONDEST_H
#define IFPACK_CONDEST_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_CondestType.h"
class Ifpack_Preconditioner;
class Epetra_RowMatrix;

double Ifpack_Condest(const Ifpack_Preconditioner& IFP,
		      const Ifpack_CondestType CT,
		      const int MaxIters = 1550,
		      const double Tol = 1e-9,
		      Epetra_RowMatrix* Matrix = 0);

#endif // IFPACK_CONDEST_H
