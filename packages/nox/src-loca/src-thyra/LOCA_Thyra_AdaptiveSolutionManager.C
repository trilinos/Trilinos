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

#include "LOCA_Thyra_AdaptiveSolutionManager.H"


using Teuchos::rcp;

LOCA::Thyra::AdaptiveSolutionManager::AdaptiveSolutionManager() :
   adaptiveMesh_(false),
   time_(0.0), iter_(0),
   p_index_(0)
{
}

void
LOCA::Thyra::AdaptiveSolutionManager::
buildSolutionGroup() {

  const NOX::Thyra::Vector initialGuess(*model_->getNominalValues().get_x());

  grp_ = Teuchos::rcp(new LOCA::Thyra::GroupWrapper(globalData_, initialGuess, model_, *paramVector_, p_index_));
  grp_->setSaveDataStrategy(saveDataStrategy_);

}

void
LOCA::Thyra::AdaptiveSolutionManager::
initialize(const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
           const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy,
           const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const Teuchos::RCP<LOCA::ParameterVector>& p,
           int p_index){

  // Create weak RCPs for the model and data strategy to address circular RCPs in the current design
  model_ = model.create_weak();
  saveDataStrategy_ = saveDataStrategy.create_weak();
  globalData_ = global_data;
  paramVector_ = p;
  p_index_ = p_index;

  buildSolutionGroup();

}

