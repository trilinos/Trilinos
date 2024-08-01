// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
