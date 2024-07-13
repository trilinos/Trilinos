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

#include "LOCA_SaveEigenData_Factory.H"
#include "LOCA_SaveEigenData_AbstractStrategy.H"
#include "LOCA_SaveEigenData_DefaultStrategy.H"

LOCA::SaveEigenData::Factory::Factory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::SaveEigenData::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>
LOCA::SaveEigenData::Factory::create(
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& eigenParams)
{
  std::string methodName = "LOCA::SaveEigenData::Factory::create()";
  Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*eigenParams);

  if (name == "Default")
    strategy =
      Teuchos::rcp(new LOCA::SaveEigenData::DefaultStrategy(globalData,
                                topParams,
                                eigenParams));
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    std::string userDefinedName = eigenParams->get(
                      "User-Defined Save Eigen Data Name",
                      "???");
    if ((*eigenParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> >(userDefinedName))
      strategy = (*eigenParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
                       methodName,
                       "Cannot find user-defined strategy: " +
                       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
                      methodName,
                      "Invalid save eigen data strategy: " +
                      name);

  return strategy;
}

const std::string&
LOCA::SaveEigenData::Factory::strategyName(
                  Teuchos::ParameterList& eigenParams) const
{
  return eigenParams.get("Save Eigen Data Method", "Default");
}
