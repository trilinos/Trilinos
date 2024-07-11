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
#include "FiniteDifference.H"
using namespace std;

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteDifference& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Vec& x, Vec& RHS)
{
  return problem.evaluate(RHS_ONLY, &x, &RHS, NULL);
}

bool Problem_Interface::computeJacobian(const Vec& x, Mat& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL, &Jac);
}

bool Problem_Interface::computePreconditioner(Mat& M)
{
  std::cout << "ERROR: Problem_Interface::computePreconditioner() - Use Explicit Jaciban only for this test problem!" << std::endl;
  throw;
}
bool Problem_Interface::preconditionVector(Vec& y)
{
  std::cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!" << std::endl;
  throw;
}
//-----------------------------------------------------------------------------

