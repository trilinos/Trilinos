// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SINGLELEVELMATLABFACTORY_DEF_HPP
#define MUELU_SINGLELEVELMATLABFACTORY_DEF_HPP
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>

#include "MueLu_Monitor.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_SingleLevelMatlabFactory_decl.hpp"
#include "MueLu_MatlabUtils_decl.hpp"

#ifdef HAVE_MUELU_MATLAB
#include "mex.h"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SingleLevelMatlabFactory()
  : hasDeclaredInput_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = getInputParamList();
  validParamList->set<std::string>("Provides", "", "A comma-separated list of objects provided by the SingleLevelMatlabFactory");
  validParamList->set<std::string>("Needs", "", "A comma-separated list of objects needed by the SingleLevelMatlabFactory");
  validParamList->set<std::string>("Function", "", "The name of the Matlab MEX function to call for Build()");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const Teuchos::ParameterList& pL = GetParameterList();
  needs_                           = tokenizeList(pL.get<std::string>("Needs"));
  // Declare inputs
  for (size_t i = 0; i < needs_.size(); i++) {
    if (!IsParamMuemexVariable(needs_[i]) && needs_[i] != "Level")
      this->Input(currentLevel, needs_[i]);
  }
  hasDeclaredInput_ = true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  const Teuchos::ParameterList& pL = GetParameterList();
  using Teuchos::rcp;
  using namespace std;
  // NOTE: mexOutput[0] is the "Provides."  Might want to modify to allow for additional outputs
  string needsList                 = pL.get<string>("Needs");
  vector<RCP<MuemexArg>> InputArgs = processNeeds<Scalar, LocalOrdinal, GlobalOrdinal, Node>(this, needsList, currentLevel);
  string providesList              = pL.get<std::string>("Provides");
  size_t numProvides               = tokenizeList(providesList).size();
  // Call mex function
  string matlabFunction = pL.get<std::string>("Function");
  if (!matlabFunction.length())
    throw std::runtime_error("Invalid matlab function name");
  vector<Teuchos::RCP<MuemexArg>> mexOutput = callMatlab(matlabFunction, numProvides, InputArgs);
  // Set output in level
  processProvides<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mexOutput, this, providesList, currentLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  const Teuchos::ParameterList& pL = GetParameterList();
  out << "SingleLevelMatlabFactory[" << pL.get<std::string>("Function") << "]";
  return out.str();
}

}  // namespace MueLu

#define MUELU_SINGLELEVELMATLABFACTORY_SHORT
#endif  // HAVE_MUELU_MATLAB

#endif  // MUELU_SINGLELEVELMATLABFACTORY_DEF_HPP
