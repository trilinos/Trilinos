// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"

#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Abstract_Factory.H"

LOCA::Factory::Factory(
      const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data),
  factory(),
  haveFactory(false),
  predictorFactory(global_data),
  continuationFactory(global_data),
  bifurcationFactory(global_data),
  stepsizeFactory(global_data),
  borderedFactory(global_data),
  eigensolverFactory(global_data),
  eigenvalueSortFactory(global_data),
  saveEigenFactory(global_data),
  anasaziOperatorFactory(global_data),
  mooreSpenceTurningPointSolverFactory(global_data),
  mooreSpencePitchforkSolverFactory(global_data),
  mooreSpenceHopfSolverFactory(global_data)
{
}

LOCA::Factory::Factory(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Abstract::Factory>& userFactory) :
  globalData(global_data),
  factory(userFactory),
  haveFactory(true),
  predictorFactory(global_data),
  continuationFactory(global_data),
  bifurcationFactory(global_data),
  stepsizeFactory(global_data),
  borderedFactory(global_data),
  eigensolverFactory(global_data),
  eigenvalueSortFactory(global_data),
  saveEigenFactory(global_data),
  anasaziOperatorFactory(global_data),
  mooreSpenceTurningPointSolverFactory(global_data),
  mooreSpencePitchforkSolverFactory(global_data),
  mooreSpenceHopfSolverFactory(global_data)
{
  // Initialize user-defined factory
  factory->init(globalData);
}

LOCA::Factory::~Factory()
{
  // Release RCP
 // if (haveFactory) factory->init(Teuchos::null);
}

Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
LOCA::Factory::createPredictorStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& predictorParams)
{
  std::string methodName = "LOCA::Factory::createPredictorStrategy()";
  Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      predictorFactory.strategyName(*predictorParams);
    bool created = factory->createPredictorStrategy(strategyName,
                            topParams,
                            predictorParams,
                            strategy);
    if (created)
      return strategy;
  }

  strategy = predictorFactory.create(topParams, predictorParams);

  return strategy;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>
LOCA::Factory::createContinuationStrategy(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& stepperParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs)
{
  std::string methodName = "LOCA::Factory::createContinuationStrategy()";
  Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      continuationFactory.strategyName(*stepperParams);
    bool created = factory->createContinuationStrategy(strategyName,
                               topParams,
                               stepperParams,
                               grp, pred, paramIDs,
                               strategy);
    if (created)
      return strategy;
  }

  strategy = continuationFactory.create(topParams, stepperParams, grp, pred,
                    paramIDs);

  return strategy;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Factory::createBifurcationStrategy(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  std::string methodName = "LOCA::Factory::createBifurcationStrategy()";
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      bifurcationFactory.strategyName(*bifurcationParams);
    bool created = factory->createBifurcationStrategy(strategyName,
                              topParams,
                              bifurcationParams,
                              grp, strategy);
    if (created)
      return strategy;
  }

  strategy = bifurcationFactory.create(topParams, bifurcationParams, grp);

  return strategy;
}

Teuchos::RCP<LOCA::StepSize::AbstractStrategy>
LOCA::Factory::createStepSizeStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams)
{
  std::string methodName = "LOCA::Factory::createStepSizeStrategy()";
  Teuchos::RCP<LOCA::StepSize::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      stepsizeFactory.strategyName(*stepsizeParams);
    bool created = factory->createStepSizeStrategy(strategyName,
                           topParams,
                           stepsizeParams,
                           strategy);
    if (created)
      return strategy;
  }

  strategy = stepsizeFactory.create(topParams, stepsizeParams);

  return strategy;
}

Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>
LOCA::Factory::createBorderedSolverStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  std::string methodName = "LOCA::Factory::createBorderedSolverStrategy()";
  Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      borderedFactory.strategyName(*solverParams);
    bool created = factory->createBorderedSolverStrategy(strategyName,
                             topParams,
                             solverParams,
                             strategy);
    if (created)
      return strategy;
  }

  strategy = borderedFactory.create(topParams, solverParams);

  return strategy;
}

Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>
LOCA::Factory::createEigensolverStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams)
{
  std::string methodName = "LOCA::Factory::createEigensolverStrategy()";
  Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      eigensolverFactory.strategyName(*eigenParams);
    bool created = factory->createEigensolverStrategy(strategyName,
                              topParams,
                              eigenParams,
                              strategy);
    if (created)
      return strategy;
  }

  strategy = eigensolverFactory.create(topParams, eigenParams);

  return strategy;
}

Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>
LOCA::Factory::createEigenvalueSortStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams)
{
  std::string methodName = "LOCA::Factory::createEigenvalueSortStrategy()";
  Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      eigenvalueSortFactory.strategyName(*eigenParams);
    bool created = factory->createEigenvalueSortStrategy(strategyName,
                             topParams,
                             eigenParams,
                             strategy);
    if (created)
      return strategy;
  }

  strategy = eigenvalueSortFactory.create(topParams, eigenParams);

  return strategy;
}

Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>
LOCA::Factory::createSaveEigenDataStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams)
{
  std::string methodName = "LOCA::Factory::createSaveEigenDataStrategy()";
  Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      saveEigenFactory.strategyName(*eigenParams);
    bool created = factory->createSaveEigenDataStrategy(strategyName,
                            topParams,
                            eigenParams,
                            strategy);
    if (created)
      return strategy;
  }

  strategy = saveEigenFactory.create(topParams, eigenParams);

  return strategy;
}

Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>
LOCA::Factory::createAnasaziOperatorStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
     const Teuchos::RCP<NOX::Abstract::Group>& grp)
{
  std::string methodName = "LOCA::Factory::createAnasaziOperatorStrategy()";
  Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      anasaziOperatorFactory.strategyName(*eigenParams);
    bool created = factory->createAnasaziOperatorStrategy(strategyName,
                              topParams,
                              eigenParams,
                              solverParams,
                              grp,
                              strategy);
    if (created)
      return strategy;
  }

  strategy = anasaziOperatorFactory.create(topParams, eigenParams,
                       solverParams, grp);

  return strategy;
}

Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy>
LOCA::Factory::createMooreSpenceTurningPointSolverStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  std::string methodName =
    "LOCA::Factory::createMooreSpenceTurningPointSolverStrategy()";
  Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      mooreSpenceTurningPointSolverFactory.strategyName(*solverParams);
    bool created =
      factory->createMooreSpenceTurningPointSolverStrategy(strategyName,
                               topParams,
                               solverParams,
                               strategy);
    if (created)
      return strategy;
  }

  strategy = mooreSpenceTurningPointSolverFactory.create(topParams,
                             solverParams);

  return strategy;
}

Teuchos::RCP<LOCA::Pitchfork::MooreSpence::SolverStrategy>
LOCA::Factory::createMooreSpencePitchforkSolverStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  std::string methodName =
    "LOCA::Factory::createMooreSpencePitchforkSolverStrategy()";
  Teuchos::RCP<LOCA::Pitchfork::MooreSpence::SolverStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      mooreSpencePitchforkSolverFactory.strategyName(*solverParams);
    bool created =
      factory->createMooreSpencePitchforkSolverStrategy(strategyName,
                            topParams,
                            solverParams,
                            strategy);
    if (created)
      return strategy;
  }

  strategy = mooreSpencePitchforkSolverFactory.create(topParams,
                              solverParams);

  return strategy;
}

Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy>
LOCA::Factory::createMooreSpenceHopfSolverStrategy(
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  std::string methodName =
    "LOCA::Factory::createMooreSpenceHopfSolverStrategy()";
  Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const std::string& strategyName =
      mooreSpenceHopfSolverFactory.strategyName(*solverParams);
    bool created =
      factory->createMooreSpenceHopfSolverStrategy(strategyName,
                           topParams,
                           solverParams,
                           strategy);
    if (created)
      return strategy;
  }

  strategy = mooreSpenceHopfSolverFactory.create(topParams,
                         solverParams);

  return strategy;
}
