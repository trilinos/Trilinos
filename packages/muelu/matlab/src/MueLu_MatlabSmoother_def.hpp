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
#include "MueLu_MatlabSmoother_decl.hpp"
#ifndef MUELU_MATLABSMOOTHER_DEF_HPP
#define MUELU_MATLABSMOOTHER_DEF_HPP
#include "MueLu_MatlabUtils_decl.hpp"

#if defined(HAVE_MUELU_MATLAB)
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MatlabSmoother(const Teuchos::ParameterList& paramList) {
  SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
  Factory::SetParameterList(paramList);
  ParameterList& pL = const_cast<ParameterList&>(this->GetParameterList());
  setupFunction_    = pL.get("Setup Function", "");
  solveFunction_    = pL.get("Solve Function", "");
  solveDataSize_    = pL.get("Number of Solver Args", 0);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  using namespace std;
  this->Input(currentLevel, "A");
  ParameterList& pL        = const_cast<ParameterList&>(this->GetParameterList());
  needsSetup_              = pL.get<string>("Needs");
  vector<string> needsList = tokenizeList(needsSetup_);
  for (size_t i = 0; i < needsList.size(); i++) {
    if (!IsParamMuemexVariable(needsList[i]) && needsList[i] != "Level")
      this->Input(currentLevel, needsList[i]);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  using namespace std;
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);
  if (this->IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::MatlabSmoother::Setup(): Setup() has already been called";
  vector<RCP<MuemexArg>> InputArgs = processNeeds<Scalar, LocalOrdinal, GlobalOrdinal, Node>(this, needsSetup_, currentLevel);
  A_                               = Factory::Get<RCP<Matrix>>(currentLevel, "A");
  RCP<MuemexArg> AmatArg           = rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<Matrix>>(A_)));
  // Always add A to the beginning of InputArgs
  InputArgs.insert(InputArgs.begin(), AmatArg);
  // Call mex function
  if (!setupFunction_.length())
    throw runtime_error("Invalid matlab function name");
  solveData_ = callMatlab(setupFunction_, solveDataSize_, InputArgs);
  this->GetOStream(Statistics1) << description() << endl;
  this->IsSetup(true);  // mark the smoother as set up
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError,
                             "MueLu::MatlabSmoother::Apply(): Setup() has not been called");
  using namespace Teuchos;
  using namespace std;
  if (InitialGuessIsZero)
    X.putScalar(0.0);
  // Push on A as first input
  vector<RCP<MuemexArg>> InputArgs;
  InputArgs.push_back(rcp(new MuemexData<RCP<Matrix>>(A_)));
  // Push on LHS & RHS
  RCP<MultiVector> Xrcp(&X, false);
  MultiVector* BPtrNonConst               = (MultiVector*)&B;
  RCP<MultiVector> Brcp                   = rcp<MultiVector>(BPtrNonConst, false);
  RCP<MuemexData<RCP<MultiVector>>> XData = rcp(new MuemexData<RCP<MultiVector>>(Xrcp));
  RCP<MuemexData<RCP<MultiVector>>> BData = rcp(new MuemexData<RCP<MultiVector>>(Brcp));
  InputArgs.push_back(XData);
  InputArgs.push_back(BData);
  for (size_t i = 0; i < solveData_.size(); i++)
    InputArgs.push_back(solveData_[i]);
  if (!solveFunction_.length()) throw std::runtime_error("Invalid matlab function name");
  vector<Teuchos::RCP<MuemexArg>> mexOutput = callMatlab(solveFunction_, 1, InputArgs);
  RCP<MuemexData<RCP<MultiVector>>> mydata  = Teuchos::rcp_static_cast<MuemexData<RCP<MultiVector>>>(mexOutput[0]);
  X                                         = *(mydata->getData());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  RCP<MatlabSmoother> smoother = rcp(new MatlabSmoother(*this));
  smoother->SetParameterList(this->GetParameterList());
  return smoother;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  if (SmootherPrototype::IsSetup()) {
    out << "Matlab Smoother(" << setupFunction_ << "/" << solveFunction_ << ")";
  } else {
    out << SmootherPrototype::description();
  }
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0)
    out << "Matlab Smoother(" << setupFunction_ << "/" << solveFunction_ << ")";

  if (verbLevel & Parameters1) {
    out0 << "Parameter list: " << std::endl;
    Teuchos::OSTab tab2(out);
    out << this->GetParameterList();
  }

  if (verbLevel & Debug) {
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }
}

// Dummy specializations for GO = long long
/*template <>
void MatlabSmoother<double,int,long long>::Setup(Level& currentLevel) {
  throw std::runtime_error("MatlabSmoother does not support GlobalOrdinal == long long.");
}
template <>
void MatlabSmoother<std::complex<double>,int,long long>::Setup(Level& currentLevel) {
  throw std::runtime_error("MatlabSmoother does not support GlobalOrdinal == long long.");
}

template <>
void MatlabSmoother<double,int,long long>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  throw std::runtime_error("MatlabSmoother does not support GlobalOrdinal == long long.");
}
template <>
void MatlabSmoother<std::complex<double>,int,long long>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  throw std::runtime_error("MatlabSmoother does not support GlobalOrdinal == long long.");
}*/

}  // namespace MueLu

#endif  // HAVE_MUELU_MATLAB
#endif  // MUELU_MATLABSMOOTHER_DEF_HPP
