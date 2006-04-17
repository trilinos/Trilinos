#include "Ifpack_ConfigDefs.h"
#include "Ifpack_DiagPreconditioner.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

// ============================================================================ 
Ifpack_DiagPreconditioner::
Ifpack_DiagPreconditioner(const Epetra_Map& DomainMap,
                          const Epetra_Map& RangeMap,
                          const Epetra_Vector& diag) : 
  UseTranspose_(false),
  DomainMap_(DomainMap),
  RangeMap_(RangeMap),
  diag_(diag)
{ }

// ============================================================================ 
Ifpack_DiagPreconditioner::~Ifpack_DiagPreconditioner()
{ }


// ============================================================================ 
int Ifpack_DiagPreconditioner::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(-1); // not defined
}

// ============================================================================ 
 int Ifpack_DiagPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-1); 

  for (int v = 0; v < X.NumVectors(); ++v)
    for (int i = 0; i < X.MyLength(); ++i)
      Y[v][i] = diag_[i] * X[v][i];
  ///Y.ReciprocalMultiply(1.0, diag_, X, 0.0);

  return(0);
}
