#ifndef IFPACK_UTILS_H
#include "Ifpack_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "unistd.h"
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_MultiVector;

//! Stops the execution of code, so that a debugger can be attached.
void Ifpack_BreakForDebugger(Epetra_Comm& Comm);

//! Creates an overlapping Epetra_CrsMatrix. Returns 0 if OverlappingLevel is 0.
Epetra_CrsMatrix* Ifpack_CreateOverlappingCrsMatrix(Epetra_RowMatrix* Matrix,
						    const int OverlappingLevel);

//! Creates an overlapping Epetra_CrsGraph. Returns 0 if OverlappingLevel is 0.
Epetra_CrsGraph* Ifpack_CreateOverlappingCrsMatrix(Epetra_CrsGraph* Graph,
						   const int OverlappingLevel);

//! Convertes an integer to string.
string Ifpack_toString(const int& x);

//! Convertes a double to string.
string Ifpack_toString(const double& x);

//! Prints on cout the true residual.
int Ifpack_PrintResidual(char* Label,  const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y);

int Ifpack_PrintResidual(const int iter, const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y);

void Ifpack_PrintSparsity(Epetra_RowMatrix& A);

#endif // IFPACK_UTILS_H
