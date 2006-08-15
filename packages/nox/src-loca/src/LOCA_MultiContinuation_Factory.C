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

#include "LOCA_MultiContinuation_Factory.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_ArcLengthGroup.H"

LOCA::MultiContinuation::Factory::Factory(
	        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::MultiContinuation::Factory::~Factory()
{
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy>
LOCA::MultiContinuation::Factory::create(
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& stepperParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs)
{
  string methodName = "LOCA::MultiContinuation::Factory::create()";
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> strategy;

  // Get name of strategy
  const string& name = strategyName(*stepperParams);

  if (name == "Natural")
    strategy = 
      Teuchos::rcp(new LOCA::MultiContinuation::NaturalGroup(globalData,
							     topParams,
							     stepperParams,
							     grp,
							     pred,
							     paramIDs));

  else if (name == "Arc Length")
    strategy = 
      Teuchos::rcp(new LOCA::MultiContinuation::ArcLengthGroup(
							      globalData,
							      topParams,
							      stepperParams,
							      grp,
							      pred,
							      paramIDs));
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = stepperParams->get("User-Defined Name",
							 "???");
    if ((*stepperParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> >(userDefinedName))
      strategy = (*stepperParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
				       methodName,
				       "Cannot find user-defined strategy: " + 
				       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
				      methodName,
				      "Invalid continuation method: " + 
				      name);

  return strategy;
}

const string&
LOCA::MultiContinuation::Factory::strategyName(
				  Teuchos::ParameterList& stepperParams) const
{
  return stepperParams.get("Continuation Method", "Arc Length");
}
