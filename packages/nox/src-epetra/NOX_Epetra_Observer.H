// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_OBSERVER
#define NOX_EPETRA_OBSERVER

#include "Epetra_Vector.h"


namespace NOX {
namespace Epetra {

class Observer
{
public:
  Observer () {};
  virtual ~Observer () {};

  //! Method called by Piro NOXSolver as a hook for postprocessing
  virtual void observeSolution(const Epetra_Vector& soln) = 0;

  //! LOCA calls this version, which will discard param info in this default implementation
  virtual void observeSolution(const Epetra_Vector& soln, double /* time_or_param_val */)
  { observeSolution(soln); };

private:

};

}
}

#endif //NOX_EPETRA_NOXOBSERVER
