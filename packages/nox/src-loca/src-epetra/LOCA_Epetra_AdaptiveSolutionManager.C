// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Epetra_AdaptiveSolutionManager.H"


using Teuchos::rcp;

LOCA::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager(
           const int number_of_time_derivatives,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_) :
   NOX::Epetra::AdaptiveSolutionManager(number_of_time_derivatives, map_, overlapMap_, overlapJacGraph_)
{
}

void
LOCA::Epetra::AdaptiveSolutionManager::
projectCurrentSolution()
{

  // Epetra_Vector element copy for now

  // Note that this must be called after the currentSolution vector is formed for the new mesh, but prior to
  // building the new solution group (grp).

  // currentSolution has been resized to hold the DOFs on the new mesh, but it is un-initialized.

  // grp->getX() is the current solution on the old mesh

  // TO provide an example, assume that the meshes are identical and we can just copy the data between them (a Copy Remesh)

    TEUCHOS_TEST_FOR_EXCEPT( currentSolution->length() != grp->getX().length());

    *currentSolution = grp->getX();

}

Teuchos::RCP<const Epetra_Vector>
LOCA::Epetra::AdaptiveSolutionManager::updateSolution(){

  // Copy new solution from group into current solution
  *currentSolution = grp->getX();

  return Teuchos::rcpFromRef(currentSolution->getEpetraVector());

}

void
LOCA::Epetra::AdaptiveSolutionManager::
getConvergenceData(int& KrylovIters, int& lastSolveKrylovIters, int& linSolves, double& tolAchieved) const {

    KrylovIters = linsys->getLinearItersTotal();
    lastSolveKrylovIters = linsys->getLinearItersLastSolve();
    linSolves = linsys->getNumLinearSolves();
    tolAchieved = linsys->getAchievedTol();

}



