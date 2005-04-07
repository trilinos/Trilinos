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
	  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	  const Teuchos::RefCountPtr<NOX::Parameter::List>& topLevelParams) :
  globalData(global_data),
  factory(),
  haveFactory(false),
  sublistParser(global_data),
  eigensolverFactory(global_data),
  eigenvalueSortFactory(global_data)
{
  reset(topLevelParams);

  // Set the factory member of the global data
  globalData->locaFactory = Teuchos::rcp(this, false);
}

LOCA::Factory::Factory(
	  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	  const Teuchos::RefCountPtr<NOX::Parameter::List>& topLevelParams,
	  const Teuchos::RefCountPtr<LOCA::Abstract::Factory>& userFactory) :
  globalData(global_data),
  factory(userFactory),
  haveFactory(true),
  sublistParser(global_data),
  eigensolverFactory(global_data),
  eigenvalueSortFactory(global_data)
{
  // Initialize user-defined factory
  factory->init(globalData);

  reset(topLevelParams);
  
  // Set the factory member of the global data
  globalData->locaFactory = Teuchos::rcp(this, false);
}

LOCA::Factory::~Factory()
{
}

NOX::Abstract::Group::ReturnType
LOCA::Factory::reset(
	    const Teuchos::RefCountPtr<NOX::Parameter::List>& topLevelParams)
{
  NOX::Abstract::Group::ReturnType result = NOX::Abstract::Group::Ok;

  // Parse sublists
  sublistParser.parseSublists(topLevelParams);

  // Reset user-provided factory if present
  if (haveFactory)
    result = factory->reset(topLevelParams);

  return result;
}

Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy>
LOCA::Factory::createEigensolverStrategy()
{
  string methodName = "LOCA::Factory::createEigensolverStrategy()";
  Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    bool created = factory->createEigensolverStrategy(strategy);
    if (created)
      return strategy;
  }

  // Get parameter lists
  Teuchos::RefCountPtr<NOX::Parameter::List> eigenParams = 
    sublistParser.getSublist("Eigensolver");
  Teuchos::RefCountPtr<NOX::Parameter::List> solverParams = 
    sublistParser.getSublist("Linear Solver");

  strategy = eigensolverFactory.create(eigenParams, solverParams);

  return strategy;
}

Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy>
LOCA::Factory::createEigenvalueSortStrategy()
{
  string methodName = "LOCA::Factory::createEigenvalueSortStrategy()";
  Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy> strategy;

  // If we have a user-provided factory, first try creating the strategy
  // using it
  if (haveFactory) {
    bool created = factory->createEigenvalueSortStrategy(strategy);
    if (created)
      return strategy;
  }

  // Get parameter lists
  Teuchos::RefCountPtr<NOX::Parameter::List> eigenParams = 
    sublistParser.getSublist("Eigensolver");

  strategy = eigenvalueSortFactory.create(eigenParams);

  return strategy;
}
