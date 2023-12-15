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
#ifndef MUELU_TWOLEVELMATLABFACTORY_DEF_HPP
#define MUELU_TWOLEVELMATLABFACTORY_DEF_HPP
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>

#include "MueLu_Monitor.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_TwoLevelMatlabFactory_decl.hpp"
#include "MueLu_MatlabUtils_decl.hpp"

#include <iostream>

#ifdef HAVE_MUELU_MATLAB
#include "mex.h"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TwoLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoLevelMatlabFactory()
  : hasDeclaredInput_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> TwoLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = getInputParamList();
  validParamList->set<std::string>("Provides", "", "A comma-separated list of objects provided on the coarse level by the TwoLevelMatlabFactory");
  validParamList->set<std::string>("Needs Fine", "", "A comma-separated list of objects needed on the fine level by the TwoLevelMatlabFactory");
  validParamList->set<std::string>("Needs Coarse", "", "A comma-separated list of objects needed on the coarse level by the TwoLevelMatlabFactory");
  validParamList->set<std::string>("Function", "", "The name of the Matlab MEX function to call for Build()");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TwoLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  using namespace std;
  const Teuchos::ParameterList& pL = GetParameterList();
  // Get needs strings
  const std::string str_nf = pL.get<std::string>("Needs Fine");
  const std::string str_nc = pL.get<std::string>("Needs Coarse");
  needsFine_               = tokenizeList(str_nf);
  needsCoarse_             = tokenizeList(str_nc);
  for (auto fineNeed : needsFine_) {
    if (!IsParamMuemexVariable(fineNeed) && fineNeed != "Level")
      this->Input(fineLevel, fineNeed);
  }
  for (auto coarseNeed : needsCoarse_) {
    if (!IsParamMuemexVariable(coarseNeed) && coarseNeed != "Level")
      this->Input(coarseLevel, coarseNeed);
  }
  hasDeclaredInput_ = true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TwoLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  const Teuchos::ParameterList& pL = GetParameterList();
  using Teuchos::rcp;
  using Teuchos::RCP;
  using namespace std;
  string needsFine                       = pL.get<string>("Needs Fine");
  string needsCoarse                     = pL.get<string>("Needs Coarse");
  vector<RCP<MuemexArg>> InputArgs       = processNeeds<Scalar, LocalOrdinal, GlobalOrdinal, Node>(this, needsFine, fineLevel);
  vector<RCP<MuemexArg>> InputArgsCoarse = processNeeds<Scalar, LocalOrdinal, GlobalOrdinal, Node>(this, needsCoarse, coarseLevel);
  // Add coarse args to the end of InputArgs
  InputArgs.reserve(InputArgs.size() + InputArgsCoarse.size());
  InputArgs.insert(InputArgs.begin(), InputArgsCoarse.begin(), InputArgsCoarse.end());

  // Determine output
  string provides    = pL.get<string>("Provides");
  size_t numProvides = tokenizeList(provides).size();
  // Call mex function
  string matlabFunction = pL.get<string>("Function");
  if (!matlabFunction.length())
    throw runtime_error("Invalid matlab function name");
  vector<RCP<MuemexArg>> mexOutput = callMatlab(matlabFunction, numProvides, InputArgs);
  processProvides<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mexOutput, this, provides, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string TwoLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  const Teuchos::ParameterList& pL = GetParameterList();
  out << "TwoLevelMatlabFactory[" << pL.get<std::string>("Function") << "]";
  return out.str();
}

}  // namespace MueLu

#define MUELU_TWOLEVELMATLABFACTORY_SHORT
#endif  // HAVE_MUELU_MATLAB

#endif  // MUELU_TWOLEVELMATLABFACTORY_DEF_HPP
