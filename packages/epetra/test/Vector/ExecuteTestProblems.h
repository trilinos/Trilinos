#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "BuildTestProblems.h"
int MatrixTests(const Epetra_BlockMap & map, const Epetra_LocalMap & LocalMap,
		    bool verbose);

int VectorTests(const Epetra_BlockMap & Map, bool verbose);

int BadResidual(bool verbose, double * Residual);
