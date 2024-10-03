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

#include "LOCA_BorderedSolver_Factory.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_BorderedSolver_Bordering.H"
#include "LOCA_BorderedSolver_Nested.H"

LOCA::BorderedSolver::Factory::Factory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::BorderedSolver::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>
LOCA::BorderedSolver::Factory::create(
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  std::string methodName = "LOCA::BorderedSolver::Factory::create()";
  Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*solverParams);

  if (name == "Bordering")
    strategy =
      Teuchos::rcp(new LOCA::BorderedSolver::Bordering(globalData,
                               topParams,
                               solverParams));

  else if (name == "Nested")
    strategy =
      Teuchos::rcp(new LOCA::BorderedSolver::Nested(globalData,
                            topParams,
                            solverParams));
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    std::string userDefinedName = solverParams->get("User-Defined Name",
                            "???");
    if ((*solverParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> >(userDefinedName))
      strategy = (*solverParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> >(userDefinedName);
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
LOCA::BorderedSolver::Factory::strategyName(
                  Teuchos::ParameterList& solverParams) const
{
  return solverParams.get("Bordered Solver Method", "Bordering");
}
