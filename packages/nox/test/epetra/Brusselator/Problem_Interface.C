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
#include "Brusselator.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(Brusselator& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec,
                       NOX::Epetra::Interface::Required::FillType fillType)
{
  return problem.evaluate(fillType, &x, &FVec, NULL);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
                    Epetra_Operator& Jac)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0, 0);
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
  std::cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << std::endl;
  throw 1;
}

bool Problem_Interface::computePreconditioner(const Epetra_Vector& x,
                          Epetra_Operator& Prec,
                          Teuchos::ParameterList* p)
{
  std::cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << std::endl;
  throw 1;
}
//-----------------------------------------------------------------------------

