// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
