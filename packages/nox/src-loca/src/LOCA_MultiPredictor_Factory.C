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
#include "LOCA_Parameter_SublistParser.H"

#include "LOCA_MultiPredictor_Factory.H"
#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "LOCA_MultiPredictor_Constant.H"
#include "LOCA_MultiPredictor_Tangent.H"
#include "LOCA_MultiPredictor_Secant.H"
#include "LOCA_MultiPredictor_Random.H"
#include "LOCA_MultiPredictor_Restart.H"

LOCA::MultiPredictor::Factory::Factory(
	        const Teuchos::RCP<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::MultiPredictor::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Factory::create(
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& predictorParams)
{
  std::string methodName = "LOCA::MultiPredictor::Factory::create()";
  Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy> strategy;

  // Get solver sublist
  Teuchos::RCP<Teuchos::ParameterList> solverParams = 
    topParams->getSublist("Linear Solver");

  // Get name of strategy
  const std::string& name = strategyName(*predictorParams);

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
    std::string userDefinedName = predictorParams->get("User-Defined Name",
							   "???");
    if ((*predictorParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy> >(userDefinedName))
      strategy = (*predictorParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy> >(userDefinedName);
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

const std::string&
LOCA::MultiPredictor::Factory::strategyName(
				  Teuchos::ParameterList& predictorParams) const
{
  return predictorParams.get("Method", "Secant");
}
