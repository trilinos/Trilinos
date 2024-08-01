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
#include "NOX_Epetra_FiniteDifferenceColoringWithUpdate.H"
#include "EpetraExt_RowMatrixOut.h"
// ----------   User Defined Includes   ----------
#include "FiniteElementProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteElementProblem& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType flag)
{
  return problem.evaluate(FiniteElementProblem::F_ONLY, &x, &FVec, NULL, flag);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  NOX::Epetra::FiniteDifferenceColoring* fdcJac = 0;
  NOX::Epetra::FiniteDifferenceColoringWithUpdate* fdcJac2 = 0;

  fdcJac = dynamic_cast<NOX::Epetra::FiniteDifferenceColoring*>(&Jac);
  if (fdcJac == 0) {
    fdcJac2 = dynamic_cast<NOX::Epetra::FiniteDifferenceColoringWithUpdate*>(&Jac);
    if(fdcJac2==0){
      std::cout << "Error: Problem_Interface::computeJacobian() - "
       << "Jacobian operator is not a NOX::Epetra::FiniteDifferenceColoring "
       << "or NOX::Epetra::FiniteDifferenceColoringWithUpdate "
       << "object!" << std::endl;
      throw "Problem_Interface Error";
    }
    else{
      fdcJac2->computeJacobian(x);
    }
  }
  else{
    fdcJac->computeJacobian(x);
  }

  return true;
}

bool Problem_Interface::computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  std::cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!" << std::endl;
  throw 1;
}
//-----------------------------------------------------------------------------

