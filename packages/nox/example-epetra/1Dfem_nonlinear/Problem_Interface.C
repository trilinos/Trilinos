// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "FiniteElementProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteElementProblem& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F)
{
  return problem.evaluate(F_ONLY, &x, &F, NULL);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  Epetra_RowMatrix* Jacobian = dynamic_cast<Epetra_RowMatrix*>(&Jac);
  if (Jacobian == NULL) {
    cout << "ERROR: Problem_Interface::computeJacobian() - The supplied"
	 << "Epetra_Operator is NOT an Epetra_RowMatrix!" << endl;
    throw;
  }
  return problem.evaluate(MATRIX_ONLY, &x, NULL, Jacobian);
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
  Epetra_RowMatrix* precMatrix = dynamic_cast<Epetra_RowMatrix*>(&M);
  if (precMatrix == NULL) {
    cout << "ERROR: Problem_Interface::computePrecMatrix() - The supplied"
	 << "Epetra_Operator is NOT an Epetra_RowMatrix!" << endl;
    throw;
  }
  return problem.evaluate(MATRIX_ONLY, &x, NULL, precMatrix);
}
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M)
{
  cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!" << endl;
  throw 1;
}
//-----------------------------------------------------------------------------

