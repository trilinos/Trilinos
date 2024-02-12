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
#ifndef MUELU_REFMAXWELLSMOOTHER_DEF_HPP
#define MUELU_REFMAXWELLSMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_RefMaxwellSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RefMaxwellSmoother(const std::string type, const Teuchos::ParameterList& paramList)
  : type_(type) {
  const bool solverSupported = (type_ == "RefMaxwell");
  this->declareConstructionOutcome(!solverSupported, "RefMaxwellSmoother does not provide the smoother '" + type_ + "'.");
  if (solverSupported)
    SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
  Factory::SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const ParameterList& pL = this->GetParameterList();

  int spaceNumber = 1;
  if (pL.isType<int>("refmaxwell: space number"))
    spaceNumber = pL.get<int>("refmaxwell: space number");
  this->Input(currentLevel, "A");
  this->Input(currentLevel, "Dk_1");
  if (spaceNumber > 1)
    this->Input(currentLevel, "Dk_2");
  this->Input(currentLevel, "D0");
  this->Input(currentLevel, "Mk_one");
  if (spaceNumber > 1)
    this->Input(currentLevel, "Mk_1_one");
  this->Input(currentLevel, "M1_beta");
  if (spaceNumber > 1)
    this->Input(currentLevel, "M1_alpha");
  this->Input(currentLevel, "invMk_1_invBeta");
  if (spaceNumber > 1)
    this->Input(currentLevel, "invMk_2_invAlpha");
  this->Input(currentLevel, "Coordinates");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);

  SetupRefMaxwell(currentLevel);
  SmootherPrototype::IsSetup(true);
  this->GetOStream(Statistics1) << description() << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupRefMaxwell(Level& currentLevel) {
  using coordinateType        = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<coordinateType, LO, GO, NO>;

  ParameterList params = this->GetParameterList();
  int spaceNumber      = 1;
  if (params.isType<int>("refmaxwell: space number"))
    spaceNumber = params.get<int>("refmaxwell: space number");

  RCP<Matrix> SM, Dk_1, Dk_2, D0, Mk_one, Mk_1_one, M1_beta, M1_alpha, invMk_1_invBeta, invMk_2_invAlpha;

  SM   = currentLevel.Get<RCP<Matrix> >("A");
  Dk_1 = currentLevel.Get<RCP<Matrix> >("Dk_1");
  if (spaceNumber > 1)
    Dk_2 = currentLevel.Get<RCP<Matrix> >("Dk_2");
  D0     = currentLevel.Get<RCP<Matrix> >("D0");
  Mk_one = currentLevel.Get<RCP<Matrix> >("Mk_one");
  if (spaceNumber > 1)
    Mk_1_one = currentLevel.Get<RCP<Matrix> >("Mk_1_one");
  M1_beta = currentLevel.Get<RCP<Matrix> >("M1_beta");
  if (spaceNumber > 1)
    M1_alpha = currentLevel.Get<RCP<Matrix> >("M1_alpha");
  invMk_1_invBeta = currentLevel.Get<RCP<Matrix> >("invMk_1_invBeta");
  if (spaceNumber > 1)
    invMk_2_invAlpha = currentLevel.Get<RCP<Matrix> >("invMk_2_invAlpha");
  auto coords = currentLevel.Get<RCP<RealValuedMultiVector> >("Coordinates");

  op_ = rcp(new MueLu::RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>(SM,
                                                                             Dk_1, Dk_2, D0,
                                                                             M1_beta, M1_alpha,
                                                                             Mk_one, Mk_1_one,
                                                                             invMk_1_invBeta, invMk_2_invAlpha,
                                                                             Teuchos::null, Teuchos::null, coords,
                                                                             params));
  std::ostringstream oss;
  op_->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(oss)));
  cachedDescription_ = oss.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::RefMaxwellSmoother::Apply(): Setup() has not been called");

  if (InitialGuessIsZero) {
    X.putScalar(0.0);
  }
  op_->apply(B, X);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  RCP<RefMaxwellSmoother> smoother = rcp(new RefMaxwellSmoother(*this));
  smoother->SetParameterList(this->GetParameterList());
  return smoother;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  if (SmootherPrototype::IsSetup()) {
    out << op_->description();
    out << std::endl
        << std::endl;
    // out << cachedDescription_;
  } else {
    out << "RefMaxwell";
  }
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters1) {
    out0 << "Parameter list: " << std::endl;
    Teuchos::OSTab tab2(out);
    out << this->GetParameterList();
  }

  if (verbLevel & External) {
    if (op_ != Teuchos::null) {
      Teuchos::OSTab tab2(out);
      out << *op_ << std::endl
          << std::endl;
    }
  }

  // if (verbLevel & Debug) {
  //   out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
  //        << "-" << std::endl
  //        << "RCP<solver_>: " << tSolver_ << std::endl;
  // }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif  // MUELU_REFMAXWELLSMOOTHER_DEF_HPP
