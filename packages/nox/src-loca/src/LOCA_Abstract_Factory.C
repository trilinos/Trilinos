// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Abstract_Factory.H"

bool
LOCA::Abstract::Factory::createPredictorStrategy(
        const std::string& /* strategyName */,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* predictorParams */,
    Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createContinuationStrategy(
    const std::string& /* strategyName */,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* stepperParams */,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& /* grp */,
    const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& /* pred */,
    const std::vector<int>& /* paramIDs */,
    Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createBifurcationStrategy(
    const std::string& /* strategyName */,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* bifurcationParams */,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& /* grp */,
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createStepSizeStrategy(
        const std::string& /* strategyName */,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* stepsizeParams */,
    Teuchos::RCP<LOCA::StepSize::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createBorderedSolverStrategy(
        const std::string& /* strategyName */,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* solverParams */,
    Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createEigensolverStrategy(
         const std::string& /* strategyName */,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
     const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */,
     Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createEigenvalueSortStrategy(
        const std::string& /* strategyName */,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */,
    Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createSaveEigenDataStrategy(
         const std::string& /* strategyName */,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
     const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */,
     Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createAnasaziOperatorStrategy(
      const std::string& /* strategyName */,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
      const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */,
      const Teuchos::RCP<Teuchos::ParameterList>& /* solverParams */,
      const Teuchos::RCP<NOX::Abstract::Group>& /* grp */,
      Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpenceTurningPointSolverStrategy(
       const std::string& /* strategyName */,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
       const Teuchos::RCP<Teuchos::ParameterList>& /* solverParams */,
       Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpencePitchforkSolverStrategy(
       const std::string& /* strategyName */,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
       const Teuchos::RCP<Teuchos::ParameterList>& /* solverParams */,
       Teuchos::RCP<LOCA::Pitchfork::MooreSpence::SolverStrategy>& /* strategy */)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpenceHopfSolverStrategy(
       const std::string& /* strategyName */,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
       const Teuchos::RCP<Teuchos::ParameterList>& /* solverParams */,
       Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy>& /* strategy */)
{
  return false;
}
