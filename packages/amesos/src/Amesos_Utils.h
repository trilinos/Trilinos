#ifndef AMESOS_UTILS_H
#include "Amesos_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "unistd.h"

class Epetra_CrsMatrix;

void Amesos_BreakForDebugger(Epetra_Comm& Comm);

Epetra_CrsMatrix* CreateOverlappingCrsMatrix(Epetra_CrsMatrix* Matrix,
					     const int OverlappingLevel);

#endif // AMESOS_UTILS_H
