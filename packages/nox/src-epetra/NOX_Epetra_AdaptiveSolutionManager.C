// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Epetra_AdaptiveSolutionManager.H"

using Teuchos::rcp;

NOX::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager(
           const int number_time_derivatives,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_) :
  num_time_deriv(number_time_derivatives)
{

   resizeMeshDataArrays(map_, overlapMap_, overlapJacGraph_);

}

void NOX::Epetra::AdaptiveSolutionManager::resizeMeshDataArrays(
           const Teuchos::RCP<const Epetra_Map> &map,
           const Teuchos::RCP<const Epetra_Map> &overlapMap,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph){

  // Create Epetra objects
  importer = rcp(new Epetra_Import(*overlapMap, *map));
  exporter = rcp(new Epetra_Export(*overlapMap, *map));
  overlapped_x = rcp(new Epetra_Vector(*overlapMap));
  if(num_time_deriv > 0)
    overlapped_xdot = rcp(new Epetra_Vector(*overlapMap));
  if(num_time_deriv > 1)
    overlapped_xdotdot = rcp(new Epetra_Vector(*overlapMap));
  overlapped_f = rcp(new Epetra_Vector(*overlapMap));
  overlapped_jac = rcp(new Epetra_CrsMatrix(Copy, *overlapJacGraph));
  tmp_ovlp_sol = rcp(new Epetra_Vector(*overlapMap));

  // Initialize solution vector and time deriv
  if(num_time_deriv > 0)
    initial_xdot = rcp(new Epetra_Vector(*map));
  if(num_time_deriv > 1)
    initial_xdotdot = rcp(new Epetra_Vector(*map));

}

//! Accessor function to Epetra_Import the solution from other PEs for output
Epetra_Vector*
NOX::Epetra::AdaptiveSolutionManager::
getOverlapSolution(const Epetra_Vector& solution) {

    tmp_ovlp_sol->Import(solution, *importer, Insert);

    return tmp_ovlp_sol.get();

}

void
NOX::Epetra::AdaptiveSolutionManager::
setInitialSolution(const Teuchos::RCP<Epetra_Vector>& soln_){

  initial_x = soln_;

  currentSolution = Teuchos::rcp(new NOX::Epetra::Vector(
    initial_x, NOX::Epetra::Vector::CreateView));

}


