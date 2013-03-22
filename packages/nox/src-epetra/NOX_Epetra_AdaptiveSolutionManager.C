//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2012) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact 
// Glen Hansen (gahanse@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Epetra_AdaptiveSolutionManager.H"

using Teuchos::rcp;

NOX::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager(
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_){


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
  overlapped_xdot = rcp(new Epetra_Vector(*overlapMap));
  overlapped_f = rcp(new Epetra_Vector(*overlapMap));
  overlapped_jac = rcp(new Epetra_CrsMatrix(Copy, *overlapJacGraph));
  tmp_ovlp_sol = rcp(new Epetra_Vector(*overlapMap));

  // Initialize solution vector and time deriv
  initial_xdot = rcp(new Epetra_Vector(*map));

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


