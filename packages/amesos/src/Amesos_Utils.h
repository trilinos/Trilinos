#ifndef AMESOS_UTILS_H
#include "Amesos_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "unistd.h"
class Epetra_CrsMatrix;

class Epetra_CrsGraph;

void Amesos_BreakForDebugger(Epetra_Comm& Comm);

Epetra_CrsMatrix* Amesos_CreateOverlappingCrsMatrix(Epetra_CrsMatrix* Matrix,
						    const int OverlappingLevel);

Epetra_CrsGraph* Amesos_CreateOverlappingCrsMatrix(Epetra_CrsGraph* Graph,
						   const int OverlappingLevel);

string Amesos_toString(const int& x);

string Amesos_toString(const double& x);

#endif // AMESOS_UTILS_H
