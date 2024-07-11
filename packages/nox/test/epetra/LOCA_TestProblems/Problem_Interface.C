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
#include "Epetra_CrsMatrix.h"
#include "Problem_Interface.H"
#include "LOCA_Parameter_Vector.H"

// ----------   User Defined Includes   ----------
#include "FiniteElementProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteElementProblem& Problem) :
  problem(Problem),
  numFillsF(0)
{
}

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag)
{
  ++numFillsF;
  return problem.evaluate(F_ONLY, &x, &F, NULL);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
                    Epetra_Operator& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,&problem.getJacobian(),
              1.0, 0.0);
}

bool Problem_Interface::computeShiftedMatrix(double alpha, double beta,
                         const Epetra_Vector& x,
                         Epetra_Operator& A)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,&problem.getJacobian(),
              alpha, beta);
}

void Problem_Interface::setParameters(const LOCA::ParameterVector& params)
{
  for (int i = 0; i < params.length(); i++ ) {
    problem.setParameter(params.getLabel(i), params.getValue(i));
  }
}
//-----------------------------------------------------------------------------

