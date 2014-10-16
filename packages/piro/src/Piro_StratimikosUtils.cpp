// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2013) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Glen Hansen (gahanse@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
  } else if (solverToken == "Rythmos") {
    if (piroParams->isSublist("Rythmos")) {
      result = Teuchos::sublist(Teuchos::sublist(piroParams, "Rythmos"), "Stratimikos");
    } else if (piroParams->isSublist("Rythmos Solver")) {
      result = Teuchos::sublist(Teuchos::sublist(piroParams, "Rythmos Solver"), "Stratimikos");
    }
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
