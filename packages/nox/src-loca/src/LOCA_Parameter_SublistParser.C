// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Parameter::SublistParser::SublistParser(
		  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) :
  globalData(global_data),
  sublistMap()
{
}

LOCA::Parameter::SublistParser::~SublistParser()
{
}

void
LOCA::Parameter::SublistParser::parseSublists(
	     const Teuchos::RefCountPtr<Teuchos::ParameterList>& topLevelParams)
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

Teuchos::RefCountPtr<Teuchos::ParameterList> 
LOCA::Parameter::SublistParser::getSublist(const string& name)
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
