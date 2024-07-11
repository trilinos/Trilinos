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
