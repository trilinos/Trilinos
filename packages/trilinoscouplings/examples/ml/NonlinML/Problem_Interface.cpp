// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
                                                                                
// ml objects
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "FiniteElementProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteElementProblem& Problem) :
problem(Problem)
{
  type_ = "Problem_Interface"; 
  return; 
}

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType flag)
{
  double t0 = GetClock();
  bool ok = problem.evaluate(F_ONLY, &x, &FVec, NULL);
  double t1 = GetClock();
  ncalls_computeF_++;
  t_ += (t1-t0);
  return ok;
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  Epetra_RowMatrix* Jacobian = dynamic_cast<Epetra_RowMatrix*>(&Jac);
  if (Jacobian == NULL) {
    std::cout << "ERROR: Problem_Interface::computeJacobian() - The supplied"
	 << "Epetra_Operator is NOT an Epetra_RowMatrix!" << std::endl;
    throw -1;
  }
  return problem.evaluate(MATRIX_ONLY, &x, NULL, Jacobian);
}

//-----------------------------------------------------------------------------
#endif // defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO)

