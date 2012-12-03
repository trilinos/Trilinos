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

#include "LOCA_Epetra_AdaptiveSolutionManager.H"


using Teuchos::rcp;

LOCA::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager(
           const Teuchos::RCP<Teuchos::ParameterList>& params,
           const Teuchos::RCP<const Epetra_Map> &map_,
           const Teuchos::RCP<const Epetra_Map> &overlapMap_,
           const Teuchos::RCP<const Epetra_CrsGraph> &overlapJacGraph_) :
   NOX::Epetra::AdaptiveSolutionManager(params, map_, overlapMap_, overlapJacGraph_)
{
}

void
LOCA::Epetra::AdaptiveSolutionManager::
projectCurrentSolution()
{

  // Epetra_Vector element copy for now

  // Note that this must be called after the currentSolution vector is formed for the new mesh, but prior to
  // building the new solution group (grp).

    *currentSolution = grp->getX();
}

void
LOCA::Epetra::AdaptiveSolutionManager::initialize(
       Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
       Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface_,
       Teuchos::RCP<Teuchos::ParameterList>& piroParams_,
       Teuchos::RCP<LOCA::ParameterVector>& pVector_,
       Teuchos::RCP<LOCA::GlobalData>& globalData_,
       bool createPrec_){

  model = model_;
  interface = interface_;
  piroParams = piroParams_;
  pVector = pVector_;
  globalData = globalData_;
  createPrec = createPrec_;

}


