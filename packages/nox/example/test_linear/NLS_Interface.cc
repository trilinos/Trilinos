// ----------   Includes   ----------
#include "NLS_Interface.h"

// ----------   User Defined Includes   ----------
#include "fill.h"

//-----------------------------------------------------------------------------
NLS_Interface::NLS_Interface()
{ }

NLS_Interface::~NLS_Interface()
{ }

void NLS_Interface::computeRHS(Epetra_Vector& x, Epetra_Vector& RHS)
{
  RHS.PutScalar(0.0);
  LO->fillMatrix(&x, &RHS, NULL);
}

void NLS_Interface::computeJacobian(Epetra_Vector& x, Epetra_RowMatrix& Jac)
{
  LO->fillMatrix(&x, NULL, &Jac);
}

void NLS_Interface::registerFill(Fill *tmp_LOPtr)
{
  LO = tmp_LOPtr;
}
//-----------------------------------------------------------------------------

