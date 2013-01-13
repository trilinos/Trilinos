// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
	        const Teuchos::RCP<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::MultiContinuation::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>
LOCA::MultiContinuation::Factory::create(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& stepperParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs)
{
  std::string methodName = "LOCA::MultiContinuation::Factory::create()";
  Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*stepperParams);

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
    std::string userDefinedName = stepperParams->get("User-Defined Name",
							 "???");
    if ((*stepperParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> >(userDefinedName))
      strategy = (*stepperParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> >(userDefinedName);
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

const std::string&
LOCA::MultiContinuation::Factory::strategyName(
				  Teuchos::ParameterList& stepperParams) const
{
  return stepperParams.get("Continuation Method", "Arc Length");
}
