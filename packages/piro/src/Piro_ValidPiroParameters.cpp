// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ValidPiroParameters.hpp"


Teuchos::RCP<const Teuchos::ParameterList>
Piro::getValidPiroParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("Valid Piro Params"));;
  validPL->sublist("NOX", false, "");
  validPL->sublist("LOCA", false, "");
  validPL->sublist("Stochastic Galerkin", false, "");
  validPL->set<std::string>("Solver Type", "","");

  // Remnant from Albany -- can be moved down a level
  validPL->sublist("VTK", false, "");

  // Maybe needs to be moved up a sublist since this 
  // isn't parsed by Piro but by the application or factory.
  validPL->set<bool>("Compute Sensitivities", false, "");

  return validPL;
}

