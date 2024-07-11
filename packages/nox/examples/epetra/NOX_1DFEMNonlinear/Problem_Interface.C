// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, const FillType flag)
{
  return problem.evaluate(F_ONLY, &x, &FVec, NULL);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  Epetra_RowMatrix* Jacobian = dynamic_cast<Epetra_RowMatrix*>(&Jac);
  if (Jacobian == NULL) {
    std::cout << "ERROR: Problem_Interface::computeJacobian() - The supplied"
     << "Epetra_Operator is NOT an Epetra_RowMatrix!" << std::endl;
    throw;
  }
  return problem.evaluate(MATRIX_ONLY, &x, NULL, Jacobian);
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
  Epetra_RowMatrix* precMatrix = dynamic_cast<Epetra_RowMatrix*>(&M);
  if (precMatrix == NULL) {
    std::cout << "ERROR: Problem_Interface::computePreconditioner() - The supplied"
     << "Epetra_Operator is NOT an Epetra_RowMatrix!" << std::endl;
    throw;
  }
  return problem.evaluate(MATRIX_ONLY, &x, NULL, precMatrix);
}
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  std::cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!" << std::endl;
  throw 1;
}
//-----------------------------------------------------------------------------

