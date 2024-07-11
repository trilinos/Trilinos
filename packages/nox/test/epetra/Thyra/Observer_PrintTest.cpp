// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Observer_PrintTest.hpp"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Solver_Generic.H"


ObserverPrintTest::ObserverPrintTest(const NOX::Utils& u) :
  numPreIterateCalls(0),
  numPostIterateCalls(0),
  numPreSolveCalls(0),
  numPostSolveCalls(0)
{
  utils = u;
}

ObserverPrintTest::~ObserverPrintTest()
{

}

void ObserverPrintTest::
runPreIterate(const NOX::Solver::Generic& solver)
{
  ++numPreIterateCalls;
  utils.out(NOX::Utils::Details) <<
    "runPreIterate() routine called!" << std::endl;
}

void ObserverPrintTest::
runPostIterate(const NOX::Solver::Generic& solver)
{
  ++numPostIterateCalls;
  utils.out(NOX::Utils::Details)
    << "runPostIterate() routine called!" << std::endl;
}

void ObserverPrintTest::
runPreSolve(const NOX::Solver::Generic& solver)
{
  ++numPreSolveCalls;
  utils.out(NOX::Utils::Details)
    << "runPreSolve() routine called!" << std::endl;
}

void ObserverPrintTest::
runPostSolve(const NOX::Solver::Generic& solver)
{
  ++numPostSolveCalls;
  utils.out(NOX::Utils::Details)
    << "runPostSolve() routine called!" << std::endl;
}
