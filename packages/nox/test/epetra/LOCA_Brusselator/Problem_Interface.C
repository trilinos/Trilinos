// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "Brusselator.H"
#include "Epetra_RowMatrix.h"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(Brusselator& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec,
                       NOX::Epetra::Interface::Required::FillType fillType)
{
  return problem.evaluate(F_ONLY, &x, &FVec, NULL);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
                    Epetra_Operator& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,
              dynamic_cast<Epetra_RowMatrix*>(&Jac),
              1.0, 0.0);
}

void Problem_Interface::setParameters(const LOCA::ParameterVector& params)
{
  problem.setParameters(params.getValue("alpha"),
            params.getValue("beta"),
            params.getValue("D1"),
            params.getValue("D2"));
}

bool Problem_Interface::computeShiftedMatrix(double alpha, double beta,
                         const Epetra_Vector& x,
                         Epetra_Operator& A)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,
              dynamic_cast<Epetra_RowMatrix*>(&A), alpha, beta);
}

//-----------------------------------------------------------------------------

