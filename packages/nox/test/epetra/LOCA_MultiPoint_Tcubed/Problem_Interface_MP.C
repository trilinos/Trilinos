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
#include "Problem_Interface_MP.H"
#include "LOCA_Parameter_Vector.H"

// ----------   User Defined Includes   ----------
#include "FiniteElementProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface_MP::Problem_Interface_MP(FiniteElementProblem& Problem) :
  problem(Problem),
  numFillsF(0)
{
}

Problem_Interface_MP::~Problem_Interface_MP()
{ }

bool Problem_Interface_MP::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag)
{
  ++numFillsF;
  return problem.evaluate(F_ONLY, &x, &F, NULL);
}

bool Problem_Interface_MP::computeJacobian(const Epetra_Vector& x,
                    Epetra_Operator& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,&problem.getJacobian());
}

void Problem_Interface_MP::setParameters(const LOCA::ParameterVector& params)
{
  for (int i = 0; i < params.length(); i++ ) {
    problem.setParameter(params.getLabel(i), params.getValue(i));
  }
}

void Problem_Interface_MP::dataForPrintSolution(const int conStep,
                          const int timeStep,
                          const int totalTimeSteps)
{
  // If this example had a printSolution method, the conStep and
  // timeStep indicies could be used to generate unique file names.
}

// The parameter that is varied in the multipoint run is set here
void Problem_Interface_MP::setMultiPointParameter(const int stepNum)
{
  problem.setParameter("Nonlinear Factor", 1.0 + 0.2*stepNum);
}

//-----------------------------------------------------------------------------

