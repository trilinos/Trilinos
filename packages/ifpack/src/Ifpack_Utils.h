#ifndef IFPACK_UTILS_H
#include "Ifpack_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "unistd.h"
class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_MultiVector;

void Ifpack_BreakForDebugger(Epetra_Comm& Comm);

Epetra_CrsMatrix* Ifpack_CreateOverlappingCrsMatrix(Epetra_CrsMatrix* Matrix,
						    const int OverlappingLevel);

Epetra_CrsGraph* Ifpack_CreateOverlappingCrsMatrix(Epetra_CrsGraph* Graph,
						   const int OverlappingLevel);

string Ifpack_toString(const int& x);

string Ifpack_toString(const double& x);

int Ifpack_PrintResidual(char* Label,  const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y);

int Ifpack_PrintResidual(int iter, const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y);

#endif // IFPACK_UTILS_H
