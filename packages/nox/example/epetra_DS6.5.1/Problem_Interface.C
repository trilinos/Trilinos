// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "DennisSchnabel.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(DennisSchnabel& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

void Problem_Interface::computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS)
{
  problem.evaluate(RHS_ONLY, &x, &RHS, NULL);
}

void Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac)
{
  problem.evaluate(MATRIX_ONLY, &x, NULL, &Jac);
}

void Problem_Interface::computePreconditioner(Epetra_RowMatrix& M)
{
  cout << "ERROR: Problem_Interface::computePreconditioner() - Use Explicit Jaciban only for this test problem!" << endl;
  throw;
}
void Problem_Interface::preconditionVector(Epetra_Vector& y)
{
  cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!" << endl;
  throw;
}
//-----------------------------------------------------------------------------

