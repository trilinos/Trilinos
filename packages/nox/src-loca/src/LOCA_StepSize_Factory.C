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
#include "LOCA_ErrorCheck.H"

#include "LOCA_StepSize_Factory.H"
#include "LOCA_StepSize_AbstractStrategy.H"
#include "LOCA_StepSize_Constant.H"
#include "LOCA_StepSize_Adaptive.H"

LOCA::StepSize::Factory::Factory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::StepSize::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::StepSize::AbstractStrategy>
LOCA::StepSize::Factory::create(
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams)
{
  std::string methodName = "LOCA::StepSize::Factory::create()";
  Teuchos::RCP<LOCA::StepSize::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*stepsizeParams);

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
    std::string userDefinedName = stepsizeParams->get("User-Defined Name",
                              "???");
    if ((*stepsizeParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::StepSize::AbstractStrategy> >(userDefinedName))
      strategy = (*stepsizeParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::StepSize::AbstractStrategy> >(userDefinedName);
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

const std::string&
LOCA::StepSize::Factory::strategyName(
                  Teuchos::ParameterList& stepsizeParams) const
{
  return stepsizeParams.get("Method", "Adaptive");
}
