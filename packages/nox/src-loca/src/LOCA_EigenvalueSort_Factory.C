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

#include "LOCA_EigenvalueSort_Factory.H"
#include "LOCA_EigenvalueSort_Strategies.H"

LOCA::EigenvalueSort::Factory::Factory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::EigenvalueSort::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>
LOCA::EigenvalueSort::Factory::create(
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams)
{
  std::string methodName = "LOCA::EigenvalueSort::Factory::create()";
  Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*eigenParams);

  if (name == "LM")
    strategy =
      Teuchos::rcp(new LOCA::EigenvalueSort::LargestMagnitude(globalData,
                                  eigenParams));

  else if (name == "LR")
    strategy =
      Teuchos::rcp(new LOCA::EigenvalueSort::LargestReal(globalData,
                             eigenParams));

  else if (name == "LI")
    strategy =
      Teuchos::rcp(new LOCA::EigenvalueSort::LargestImaginary(globalData,
                                  eigenParams));

  else if (name == "SM")
    strategy =
      Teuchos::rcp(new LOCA::EigenvalueSort::SmallestMagnitude(globalData,
                                   eigenParams));

  else if (name == "SR")
    strategy =
      Teuchos::rcp(new LOCA::EigenvalueSort::SmallestReal(globalData,
                              eigenParams));

  else if (name == "SI")
    strategy =
      Teuchos::rcp(new LOCA::EigenvalueSort::SmallestImaginary(globalData,
                                   eigenParams));

  else if (name == "CA")
    strategy =
      Teuchos::rcp(new LOCA::EigenvalueSort::LargestRealInverseCayley(
                                   globalData,
                                   eigenParams));

  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    std::string userDefinedName =
      eigenParams->get("User-Defined Sorting Method Name",
                "???");
    if ((*eigenParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy> >(userDefinedName))
      strategy = (*eigenParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
                 methodName,
                 "Cannot find user-defined sorting strategy: " +
                 userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
                      methodName,
                      "Invalid sorting strategy: " +
                      name);

  return strategy;
}

const std::string&
LOCA::EigenvalueSort::Factory::strategyName(
                  Teuchos::ParameterList& eigenParams) const
{
  return eigenParams.get("Sorting Order", "LM");
}
