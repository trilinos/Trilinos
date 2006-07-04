// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "Teuchos_ParameterList.hpp"

#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Abstract_Factory.H"

LOCA::Factory::Factory(
	  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) :
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
  mooreSpencePitchforkSolverFactory(global_data)
{
  // Set the factory member of the global data
  globalData->locaFactory = Teuchos::rcp(this, false);
}

LOCA::Factory::Factory(
	  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	  const Teuchos::RefCountPtr<LOCA::Abstract::Factory>& userFactory) :
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
  mooreSpencePitchforkSolverFactory(global_data)
{
  // Initialize user-defined factory
  factory->init(globalData);
  
  // Set the factory member of the global data
  globalData->locaFactory = Teuchos::rcp(this, false);
}

LOCA::Factory::~Factory()
{
}

Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
LOCA::Factory::createPredictorStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& predictorParams)
{
  string methodName = "LOCA::Factory::createPredictorStrategy()";
  Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy>
LOCA::Factory::createContinuationStrategy(
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& stepperParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs)
{
  string methodName = "LOCA::Factory::createContinuationStrategy()";
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::Factory::createBifurcationStrategy(
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  string methodName = "LOCA::Factory::createBifurcationStrategy()";
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy>
LOCA::Factory::createStepSizeStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& stepsizeParams)
{
  string methodName = "LOCA::Factory::createStepSizeStrategy()";
  Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy>
LOCA::Factory::createBorderedSolverStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams)
{
  string methodName = "LOCA::Factory::createBorderedSolverStrategy()";
  Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy>
LOCA::Factory::createEigensolverStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams)
{
  string methodName = "LOCA::Factory::createEigensolverStrategy()";
  Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy>
LOCA::Factory::createEigenvalueSortStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams)
{
  string methodName = "LOCA::Factory::createEigenvalueSortStrategy()";
  Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::SaveEigenData::AbstractStrategy>
LOCA::Factory::createSaveEigenDataStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams)
{
  string methodName = "LOCA::Factory::createSaveEigenDataStrategy()";
  Teuchos::RefCountPtr<LOCA::SaveEigenData::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy>
LOCA::Factory::createAnasaziOperatorStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
	 const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp)
{
  string methodName = "LOCA::Factory::createAnasaziOperatorStrategy()";
  Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::SolverStrategy>
LOCA::Factory::createMooreSpenceTurningPointSolverStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams)
{
  string methodName = 
    "LOCA::Factory::createMooreSpenceTurningPointSolverStrategy()";
  Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::SolverStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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

Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::SolverStrategy>
LOCA::Factory::createMooreSpencePitchforkSolverStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams)
{
  string methodName = 
    "LOCA::Factory::createMooreSpencePitchforkSolverStrategy()";
  Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::SolverStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    const string& strategyName = 
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
