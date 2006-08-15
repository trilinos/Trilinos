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

#include "LOCA_StepSize_Factory.H"
#include "LOCA_StepSize_AbstractStrategy.H"
#include "LOCA_StepSize_Constant.H"
#include "LOCA_StepSize_Adaptive.H"

LOCA::StepSize::Factory::Factory(
	        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::StepSize::Factory::~Factory()
{
}

Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy>
LOCA::StepSize::Factory::create(
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& stepsizeParams)
{
  string methodName = "LOCA::StepSize::Factory::create()";
  Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy> strategy;

  // Get name of strategy
  const string& name = strategyName(*stepsizeParams);

  if (name == "Constant")
    strategy = 
      Teuchos::rcp(new LOCA::StepSize::Constant(globalData,
						topParams,
						stepsizeParams));
  else if (name == "Adaptive")
    strategy = 
      Teuchos::rcp(new LOCA::StepSize::Adaptive(globalData,
						topParams,
						stepsizeParams));
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = stepsizeParams->get("User-Defined Name",
							  "???");
    if ((*stepsizeParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy> >(userDefinedName))
      strategy = (*stepsizeParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
				       methodName,
				       "Cannot find user-defined strategy: " + 
				       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
				      methodName,
				      "Invalid step size control strategy: " + 
				      name);

  return strategy;
}

const string&
LOCA::StepSize::Factory::strategyName(
				  Teuchos::ParameterList& stepsizeParams) const
{
  return stepsizeParams.get("Method", "Adaptive");
}
