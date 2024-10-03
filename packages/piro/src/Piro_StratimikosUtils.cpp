// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_StratimikosUtils.hpp"

Teuchos::RCP<Teuchos::ParameterList>
Piro::extractStratimikosParams(const Teuchos::RCP<Teuchos::ParameterList> &piroParams)
{
  Teuchos::RCP<Teuchos::ParameterList> result;

  const std::string solverToken = piroParams->get<std::string>("Solver Type");
  if (solverToken == "NOX" || solverToken == "LOCA" || solverToken == "LOCA Adaptive") {
    result = Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(
                piroParams, "NOX"), "Direction"), "Newton"), "Stratimikos Linear Solver"), "Stratimikos");
  } 
  else if (solverToken == "Tempus") {
    result = Teuchos::sublist(Teuchos::sublist(piroParams, "Tempus"), "Stratimikos");
  } 
  else if (solverToken == "Trapezoid Rule") {
    result = Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(
                piroParams, "NOX"), "Direction"), "Newton"), "Stratimikos Linear Solver"), "Stratimikos");
  }

  return result;
}

void
Piro::renamePreconditionerParamList(const Teuchos::RCP<Teuchos::ParameterList> &stratParams, 
      const std::string &oldname, const std::string &newname){


  if(stratParams->getPtr<std::string>("Preconditioner Type")){

     const std::string currentval = stratParams->get<std::string>("Preconditioner Type");

     if(currentval == oldname)

       stratParams->set<std::string>("Preconditioner Type", newname);

     else

       return; // do nothing if the names do not match

  }
  else 
     return; // return if Preconditioner Type isn't specified

  // Does the old sublist exist?
  if (stratParams->isSublist("Preconditioner Types") && stratParams->sublist("Preconditioner Types").isSublist(oldname)) {
      Teuchos::ParameterList &result = stratParams->sublist("Preconditioner Types").sublist(oldname, true);
      result.setName(newname);
  }

}
