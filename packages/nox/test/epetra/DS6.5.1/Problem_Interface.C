// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ----------   Includes   ----------
#include "Problem_Interface.H"
#include "DennisSchnabel.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(DennisSchnabel& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec,
                     NOX::Epetra::Interface::Required::FillType fillType)
{
  return problem.evaluate(fillType, &x, &FVec);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
                    Epetra_Operator& Jac)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, NULL);
}

//-----------------------------------------------------------------------------

