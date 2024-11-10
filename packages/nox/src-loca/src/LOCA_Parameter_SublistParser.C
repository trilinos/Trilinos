// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Parameter::SublistParser::SublistParser(
          const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data),
  sublistMap()
{
}

LOCA::Parameter::SublistParser::~SublistParser()
{
}

void
LOCA::Parameter::SublistParser::parseSublists(
         const Teuchos::RCP<Teuchos::ParameterList>& topLevelParams)
{
  // Top level sublist
  sublistMap["Top Level"] = topLevelParams;

  // LOCA sublist
  Teuchos::ParameterList& locaSublist = topLevelParams->sublist("LOCA");
  sublistMap["LOCA"] = Teuchos::rcp(&locaSublist, false);

  // Stepper sublist
  Teuchos::ParameterList& stepperSublist = locaSublist.sublist("Stepper");
  sublistMap["Stepper"] = Teuchos::rcp(&stepperSublist, false);

  // Eigensolver sublist
  Teuchos::ParameterList& eigensolverSublist =
    stepperSublist.sublist("Eigensolver");
  sublistMap["Eigensolver"] = Teuchos::rcp(&eigensolverSublist, false);

  // Constraints sublist
  Teuchos::ParameterList& constraintsSublist =
    locaSublist.sublist("Constraints");
  sublistMap["Constraints"] = Teuchos::rcp(&constraintsSublist, false);

  // Bifurcation sublist
  Teuchos::ParameterList& bifurcationSublist =
    locaSublist.sublist("Bifurcation");
  sublistMap["Bifurcation"] = Teuchos::rcp(&bifurcationSublist, false);

  // Predictor sublist
  Teuchos::ParameterList& predictorSublist = locaSublist.sublist("Predictor");
  sublistMap["Predictor"] = Teuchos::rcp(&predictorSublist, false);

  // First Step Predictor sublist
  Teuchos::ParameterList& fspredictorSublist =
    predictorSublist.sublist("First Step Predictor");
  sublistMap["First Step Predictor"] =
    Teuchos::rcp(&fspredictorSublist, false);

  // Last Step Predictor sublist
  Teuchos::ParameterList& lspredictorSublist =
    predictorSublist.sublist("Last Step Predictor");
  sublistMap["Last Step Predictor"] = Teuchos::rcp(&lspredictorSublist, false);

  // Stepsize sublist
  Teuchos::ParameterList& stepsizeSublist = locaSublist.sublist("Step Size");
  sublistMap["Step Size"] = Teuchos::rcp(&stepsizeSublist, false);

  // NOX sublist
  Teuchos::ParameterList& noxSublist = topLevelParams->sublist("NOX");
  sublistMap["NOX"] = Teuchos::rcp(&noxSublist, false);

  // Direction sublist
  Teuchos::ParameterList& directionSublist = noxSublist.sublist("Direction");
  sublistMap["Direction"] = Teuchos::rcp(&directionSublist, false);

  // Newton sublist
  Teuchos::ParameterList& newtonSublist = directionSublist.sublist("Newton");
  sublistMap["Newton"] = Teuchos::rcp(&newtonSublist, false);

  // Linear Solver sublist
  Teuchos::ParameterList& lsSublist = newtonSublist.sublist("Linear Solver");
  sublistMap["Linear Solver"] = Teuchos::rcp(&lsSublist, false);

  // Line Search sublist
  Teuchos::ParameterList& lineSearchSublist = noxSublist.sublist("Line Search");
  sublistMap["Line Search"] = Teuchos::rcp(&lineSearchSublist, false);

  // Printing sublist
  Teuchos::ParameterList& printingSublist = noxSublist.sublist("Printing");
  sublistMap["Printing"] = Teuchos::rcp(&printingSublist, false);
}

Teuchos::RCP<Teuchos::ParameterList>
LOCA::Parameter::SublistParser::getSublist(const std::string& name)
{
  // Find name in list, if it exists.
  SublistMapIterator i = sublistMap.find(name);

  // If it does not exist throw an error.
  if (i == sublistMap.end()) {
   globalData->locaErrorCheck->throwError(
                 "LOCA::Parameter::SublistParser::getSublist()",
                 "Invalid sublist name: " + name);
  }

  // Return sublist
  return (*i).second;
}
