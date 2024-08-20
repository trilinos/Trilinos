// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_ADAPTMANAGER
#define NOX_EPETRA_ADAPTMANAGER

#include "Epetra_Vector.h"


namespace NOX {
namespace Epetra {

class AdaptManager
{
public:

  AdaptManager () :
    time(0.0), iter(0) {};

  virtual ~AdaptManager () {};

  //! Method called by Piro NOXSolver to determine if the mesh needs adapting
  // A return type of true means that the mesh should be adapted
  virtual bool queryAdaptationCriteria() = 0;

  //! Apply adaptation method to mesh and problem. Returns true if adaptation is performed successfully.
  virtual bool adaptMesh(const Epetra_Vector& solution, const Epetra_Vector& ovlp_solution) = 0;

  //! Transfer solution between meshes.
  virtual void solutionTransfer(const Epetra_Vector& oldSolution,
        Epetra_Vector& newSolution) = 0;

  //! Track the time state of the mesh
  void setTime(const double time_){ time = time_; }

  //! Track the time state of the mesh
  void setIteration(const int iter_){ iter = iter_; }

protected:

  double time;
  int iter;

};

}
}

#endif //NOX_EPETRA_NOXADAPTMANAGER
