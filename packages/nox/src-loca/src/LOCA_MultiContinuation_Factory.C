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
