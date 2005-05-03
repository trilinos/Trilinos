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

#include "NOX_Parameter_List.H"

#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Abstract_Factory.H"

// Factories
#include "LOCA_Eigensolver_Factory.H"

LOCA::Factory::Factory(
	  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) :
  globalData(global_data),
  factory(),
  haveFactory(false),
  predictorFactory(global_data),
  continuationFactory(global_data),
  borderedFactory(global_data),
  eigensolverFactory(global_data),
  eigenvalueSortFactory(global_data)
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
  borderedFactory(global_data),
  eigensolverFactory(global_data),
  eigenvalueSortFactory(global_data)
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
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& predictorParams)
{
  string methodName = "LOCA::Factory::createPredictorStrategy()";
  Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    bool created = factory->createPredictorStrategy(topParams,
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
      const Teuchos::RefCountPtr<NOX::Parameter::List>& stepperParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs)
{
  string methodName = "LOCA::Factory::createContinuationStrategy()";
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    bool created = factory->createContinuationStrategy(topParams,
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

Teuchos::RefCountPtr<LOCA::BorderedSystem::AbstractStrategy>
LOCA::Factory::createBorderedSystemStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams)
{
  string methodName = "LOCA::Factory::createBorderedSystemStrategy()";
  Teuchos::RefCountPtr<LOCA::BorderedSystem::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    bool created = factory->createBorderedSystemStrategy(topParams,
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
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams)
{
  string methodName = "LOCA::Factory::createEigensolverStrategy()";
  Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    bool created = factory->createEigensolverStrategy(topParams,
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
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams)
{
  string methodName = "LOCA::Factory::createEigenvalueSortStrategy()";
  Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    bool created = factory->createEigenvalueSortStrategy(topParams,
							 eigenParams,
							 strategy);
    if (created)
      return strategy;
  }

  strategy = eigenvalueSortFactory.create(topParams, eigenParams);

  return strategy;
}
