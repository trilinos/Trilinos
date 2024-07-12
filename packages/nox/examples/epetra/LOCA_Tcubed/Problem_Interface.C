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
  paramLib()
{
  paramLib.addParameterEntry("Nonlinear Factor", problem,
                 &FiniteElementProblem::factor);
  paramLib.addParameterEntry("Left BC", problem,
                 &FiniteElementProblem::leftBC);
  paramLib.addParameterEntry("Right BC", problem,
                 &FiniteElementProblem::rightBC);
}

Problem_Interface::~Problem_Interface()
{
}

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType flag)
{
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
    //problem.set(params.getLabel(i), params.getValue(i));
    paramLib.setValue(params.getLabel(i), params.getValue(i));
  }
}

void Problem_Interface::printSolution(const Epetra_Vector& x, double conParam)
{
  double ave; x.MeanValue(&ave);
  std::cout << "   FreeEnergy Mock Function (param, soln norm, FreeEnergy)= "
       << conParam << "   " <<  ave << "   " << computeFreeEnergy(x) << std::endl;

  problem.printSolution(x, conParam);
}

double Problem_Interface::computeFreeEnergy(const Epetra_Vector& x)
{
    double ave; x.MeanValue(&ave);
    return abs(ave - 1.2);
}
//-----------------------------------------------------------------------------

