// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"

LOCA::GlobalData::GlobalData(
           const Teuchos::RCP<NOX::Utils>& loca_utils,
           const Teuchos::RCP<LOCA::ErrorCheck>& loca_error_check,
           const Teuchos::RCP<LOCA::Factory>& loca_factory) :
  locaUtils(loca_utils),
  locaErrorCheck(loca_error_check),
  locaFactory(loca_factory),
  parsedParams()
{
}

LOCA::GlobalData::~GlobalData()
{
}

Teuchos::RCP<LOCA::GlobalData>
LOCA::createGlobalData(
          const Teuchos::RCP<Teuchos::ParameterList>& paramList,
          const Teuchos::RCP<LOCA::Abstract::Factory>& userFactory)
{
  // Create a global data object with null data fields
  Teuchos::RCP<LOCA::GlobalData> globalData =
    Teuchos::rcp(new LOCA::GlobalData(Teuchos::null,
                      Teuchos::null,
                      Teuchos::null));

  // Create utils
  globalData->locaUtils =
    Teuchos::rcp(new NOX::Utils(paramList->sublist("NOX").sublist("Printing")));

  // Create error check
  globalData->locaErrorCheck =
    Teuchos::rcp(new LOCA::ErrorCheck(globalData));

  // Create factory
  if (userFactory != Teuchos::null)
    globalData->locaFactory = Teuchos::rcp(new LOCA::Factory(globalData,
                                 userFactory));
  else
    globalData->locaFactory = Teuchos::rcp(new LOCA::Factory(globalData));

  // Parse parameter list
  globalData->parsedParams =
    Teuchos::rcp(new Parameter::SublistParser(globalData));
  globalData->parsedParams->parseSublists(paramList);

  return globalData;
}

void
LOCA::destroyGlobalData(
            const Teuchos::RCP<LOCA::GlobalData>& globalData)
{
  globalData->locaUtils = Teuchos::null;
  globalData->locaErrorCheck = Teuchos::null;
  globalData->locaFactory = Teuchos::null;
  globalData->parsedParams = Teuchos::null;
}
