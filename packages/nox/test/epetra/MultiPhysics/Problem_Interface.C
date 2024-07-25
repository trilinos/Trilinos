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
#include "GenericEpetraProblem.H"
#include "Problem_Manager.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(GenericEpetraProblem& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType flag)
{
  return problem.evaluate(flag, &x, &FVec);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
                    Epetra_Operator& Jac)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0);
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Prec, &x, 0);
}
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x,
                          Epetra_Operator& Prec,
                       Teuchos::ParameterList* precParams)
{
  // Pass through to let the Problem fill its owned matrix
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0);
}
//-----------------------------------------------------------------------------

