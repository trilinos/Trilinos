// ----------   Includes   ----------
#include "Problem_Interface.h"

// ----------   User Defined Includes   ----------
#include "fill.h"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface()
{ }

Problem_Interface::~Problem_Interface()
{ }

void Problem_Interface::computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS)
{
  RHS.PutScalar(0.0);
  LO->fillMatrix(&x, &RHS, NULL);
}

void Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac)
{
  LO->fillMatrix(&x, NULL, &Jac);
}

void Problem_Interface::registerFill(Fill *tmp_LOPtr)
{
  LO = tmp_LOPtr;
}
//-----------------------------------------------------------------------------

