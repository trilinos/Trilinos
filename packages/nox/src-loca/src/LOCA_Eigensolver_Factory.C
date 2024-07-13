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

#include "LOCA_Eigensolver_Factory.H"
#include "LOCA_Eigensolver_AbstractStrategy.H"
#include "LOCA_Eigensolver_DefaultStrategy.H"
#ifdef HAVE_LOCA_ANASAZI
#include "LOCA_Eigensolver_AnasaziStrategy.H"
#endif

LOCA::Eigensolver::Factory::Factory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::Eigensolver::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>
LOCA::Eigensolver::Factory::create(
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& eigenParams)
{
  std::string methodName = "LOCA::Eigensolver::Factory::create()";
  Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*eigenParams);

  if (name == "Default")
    strategy =
      Teuchos::rcp(new LOCA::Eigensolver::DefaultStrategy(globalData,
                              topParams,
                              eigenParams));
  else if (name == "Anasazi") {
#ifdef HAVE_LOCA_ANASAZI
    strategy =
      Teuchos::rcp(new LOCA::Eigensolver::AnasaziStrategy(globalData,
                              topParams,
                              eigenParams));
#else
    globalData->locaErrorCheck->throwError(methodName,
                       "Anasazi strategy requested, but LOCA was not configured with Anasazi support enabled.");
#endif
  }
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    std::string userDefinedName = eigenParams->get("User-Defined Name",
                               "???");
    if ((*eigenParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy> >(userDefinedName))
      strategy = (*eigenParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
                       methodName,
                       "Cannot find user-defined strategy: " +
                       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
                      methodName,
                      "Invalid eigensolver strategy: " +
                      name);

  return strategy;
}

const std::string&
LOCA::Eigensolver::Factory::strategyName(
                  Teuchos::ParameterList& eigenParams) const
{
  return eigenParams.get("Method", "Default");
}
