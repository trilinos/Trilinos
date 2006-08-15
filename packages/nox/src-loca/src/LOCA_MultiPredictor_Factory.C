// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Parameter_SublistParser.H"

#include "LOCA_MultiPredictor_Factory.H"
#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "LOCA_MultiPredictor_Constant.H"
#include "LOCA_MultiPredictor_Tangent.H"
#include "LOCA_MultiPredictor_Secant.H"
#include "LOCA_MultiPredictor_Random.H"
#include "LOCA_MultiPredictor_Restart.H"

LOCA::MultiPredictor::Factory::Factory(
	        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::MultiPredictor::Factory::~Factory()
{
}

Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Factory::create(
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& predictorParams)
{
  string methodName = "LOCA::MultiPredictor::Factory::create()";
  Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> strategy;

  // Get solver sublist
  Teuchos::RefCountPtr<Teuchos::ParameterList> solverParams = 
    topParams->getSublist("Linear Solver");

  // Get name of strategy
  const string& name = strategyName(*predictorParams);

  if (name == "Constant")
    strategy = 
      Teuchos::rcp(new LOCA::MultiPredictor::Constant(globalData,
						      predictorParams));

  else if (name == "Tangent")
    strategy = 
      Teuchos::rcp(new LOCA::MultiPredictor::Tangent(globalData,
						     predictorParams,
						     solverParams));
  else if (name == "Secant")
    strategy = 
      Teuchos::rcp(new LOCA::MultiPredictor::Secant(globalData,
						    topParams,
						    predictorParams));
  else if (name == "Random")
    strategy = 
      Teuchos::rcp(new LOCA::MultiPredictor::Random(globalData,
						    predictorParams));
  else if (name == "Restart")
    strategy = 
      Teuchos::rcp(new LOCA::MultiPredictor::Restart(globalData,
						    predictorParams));

  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = predictorParams->get("User-Defined Name",
							   "???");
    if ((*predictorParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> >(userDefinedName))
      strategy = (*predictorParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
				       methodName,
				       "Cannot find user-defined strategy: " + 
				       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
				      methodName,
				      "Invalid predictor strategy: " + 
				      name);

  return strategy;
}

const string&
LOCA::MultiPredictor::Factory::strategyName(
				  Teuchos::ParameterList& predictorParams) const
{
  return predictorParams.get("Method", "Secant");
}
