// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_AdaptiveSolutionManager.hpp"


using Teuchos::rcp;

Thyra::AdaptiveSolutionManager::AdaptiveSolutionManager() :
   adaptiveMesh_(false),
   time_(0.0), iter_(0)
{
}

void
Thyra::LOCAAdaptiveState::
buildSolutionGroup() {

  const NOX::Thyra::Vector initialGuess(*model_->getNominalValues().get_x());

  grp_ = Teuchos::rcp(new LOCA::Thyra::GroupWrapper(globalData_, initialGuess, model_, *paramVector_, p_index_));
  grp_->setSaveDataStrategy(saveDataStrategy_);

}

Thyra::LOCAAdaptiveState::
LOCAAdaptiveState(const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
           const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy,
           const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const Teuchos::RCP<LOCA::ParameterVector>& p,
           int p_index)
    : AdaptiveStateBase(model),
      saveDataStrategy_(saveDataStrategy.create_weak()),
      globalData_(global_data),
      paramVector_(p),
      p_index_(p_index)
{

  // Create weak RCPs for the model and data strategy to address circular RCPs in the current design

  buildSolutionGroup();

}

Thyra::TransAdaptiveState::
TransAdaptiveState(const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model)
    : AdaptiveStateBase(model)
{
}

Thyra::AdaptiveStateBase::
AdaptiveStateBase(const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model)
   : model_(model.create_weak())
{

  // Create weak RCPs for the model and data strategy to address circular RCPs in the current design

}

