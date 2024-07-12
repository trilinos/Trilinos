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

#include "LOCA_TurningPoint_MooreSpence_SolverFactory.H"
#include "LOCA_TurningPoint_MooreSpence_SolverStrategy.H"
#include "LOCA_TurningPoint_MooreSpence_SalingerBordering.H"
#include "LOCA_TurningPoint_MooreSpence_PhippsBordering.H"

LOCA::TurningPoint::MooreSpence::SolverFactory::SolverFactory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::TurningPoint::MooreSpence::SolverFactory::~SolverFactory()
{
}

Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy>
LOCA::TurningPoint::MooreSpence::SolverFactory::create(
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  std::string methodName =
    "LOCA::TurningPoint::MooreSpence::SolverFactory::create()";
  Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*solverParams);

  if (name == "Salinger Bordering")
    strategy =
      Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::SalingerBordering(
                                globalData,
                                topParams,
                                solverParams));

  else if (name == "Phipps Bordering")
    strategy =
      Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::PhippsBordering(
                                globalData,
                                topParams,
                                solverParams));

  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    std::string userDefinedName = solverParams->get("User-Defined Name",
                            "???");
    if ((*solverParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy> >(userDefinedName))
      strategy = (*solverParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
                       methodName,
                       "Cannot find user-defined strategy: " +
                       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
                      methodName,
                      "Invalid bordered solver strategy: " +
                      name);

  return strategy;
}

const std::string&
LOCA::TurningPoint::MooreSpence::SolverFactory::strategyName(
                  Teuchos::ParameterList& solverParams) const
{
  return solverParams.get("Solver Method", "Salinger Bordering");
}
