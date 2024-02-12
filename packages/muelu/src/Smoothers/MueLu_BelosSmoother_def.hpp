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
#ifndef MUELU_BELOSSMOOTHER_DEF_HPP
#define MUELU_BELOSSMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_BELOS)

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMultiVector.hpp>
#endif

#include "MueLu_BelosSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BelosSmoother(const std::string type, const Teuchos::ParameterList& paramList)
  : type_(type) {
  bool solverSupported = false;
  Belos::SolverFactory<Scalar, tMV, tOP> tFactory;
  solverSupported = solverSupported || tFactory.isSupported(type);
  // if (!solverSupported) {
  //   Teuchos::Array<std::string> supportedSolverNames = factory.supportedSolverNames();

  //   std::ostringstream outString;
  //   outString << "[";
  //   for (auto iter = supportedSolverNames.begin(); iter != supportedSolverNames.end(); ++iter) {
  //     outString << "\"" << *iter << "\"";
  //     if (iter + 1 != supportedSolverNames.end()) {
  //       outString << ", ";
  //     }
  //   }
  //   outString << "]";

  //   TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Belos does not provide the solver '" << type_ << "'.\nSupported Solvers: " << outString.str());
  // }
  this->declareConstructionOutcome(!solverSupported, "Belos does not provide the smoother '" + type_ + "'.");
  if (solverSupported)
    SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
  Factory::SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  this->Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);

  A_ = Factory::Get<RCP<Matrix> >(currentLevel, "A");
  SetupBelos(currentLevel);
  SmootherPrototype::IsSetup(true);
  this->GetOStream(Statistics1) << description() << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupBelos(Level& currentLevel) {
  bool useTpetra = A_->getRowMap()->lib() == Xpetra::UseTpetra;

  if (useTpetra) {
    tBelosProblem_ = rcp(new Belos::LinearProblem<Scalar, tMV, tOP>());
    RCP<tOP> tA    = Utilities::Op2NonConstTpetraCrs(A_);
    tBelosProblem_->setOperator(tA);

    Belos::SolverFactory<SC, tMV, tOP> solverFactory;
    tSolver_ = solverFactory.create(type_, rcpFromRef(const_cast<ParameterList&>(this->GetParameterList())));
    tSolver_->setProblem(tBelosProblem_);
  } else {
    TEUCHOS_ASSERT(false);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::BelosSmoother::Apply(): Setup() has not been called");

  if (A_->getRowMap()->lib() == Xpetra::UseTpetra) {
    if (InitialGuessIsZero) {
      X.putScalar(0.0);

      RCP<Tpetra::MultiVector<SC, LO, GO, NO> > tpX       = rcpFromRef(Utilities::MV2NonConstTpetraMV(X));
      RCP<const Tpetra::MultiVector<SC, LO, GO, NO> > tpB = rcpFromRef(Utilities::MV2TpetraMV(B));

      tBelosProblem_->setInitResVec(tpB);
      tBelosProblem_->setProblem(tpX, tpB);
      tSolver_->solve();

    } else {
      typedef Teuchos::ScalarTraits<Scalar> TST;
      RCP<MultiVector> Residual   = Utilities::Residual(*A_, X, B);
      RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());

      RCP<Tpetra::MultiVector<SC, LO, GO, NO> > tpX       = rcpFromRef(Utilities::MV2NonConstTpetraMV(*Correction));
      RCP<const Tpetra::MultiVector<SC, LO, GO, NO> > tpB = rcpFromRef(Utilities::MV2TpetraMV(*Residual));

      tBelosProblem_->setInitResVec(tpB);
      tBelosProblem_->setProblem(tpX, tpB);
      tSolver_->solve();

      X.update(TST::one(), *Correction, TST::one());
    }
  } else {
    TEUCHOS_ASSERT(false);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  RCP<BelosSmoother> smoother = rcp(new BelosSmoother(*this));
  smoother->SetParameterList(this->GetParameterList());
  return smoother;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  if (SmootherPrototype::IsSetup()) {
    if (A_->getRowMap()->lib() == Xpetra::UseTpetra) {
      out << tSolver_->description();
    }
  } else {
    out << "BELOS {type = " << type_ << "}";
  }
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters1) {
    out0 << "Parameter list: " << std::endl;
    Teuchos::OSTab tab2(out);
    out << this->GetParameterList();
  }

  if (verbLevel & External) {
    if (tSolver_ != Teuchos::null) {
      Teuchos::OSTab tab2(out);
      out << *tSolver_ << std::endl;
    }
  }

  if (verbLevel & Debug) {
    if (A_->getRowMap()->lib() == Xpetra::UseTpetra) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<solver_>: " << tSolver_ << std::endl;
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif  // HAVE_MUELU_BELOS
#endif  // MUELU_BELOSSMOOTHER_DEF_HPP
