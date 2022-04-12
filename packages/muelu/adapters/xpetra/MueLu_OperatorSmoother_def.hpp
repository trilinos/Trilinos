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
#ifndef MUELU_OperatorSMOOTHER_DEF_HPP
#define MUELU_OperatorSMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_OperatorSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"



namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::OperatorSmoother(const std::string type, const Teuchos::ParameterList& paramList)
  : type_(type)
  {
    const bool solverSupported = (type_ == "RefMaxwell");
    this->declareConstructionOutcome(!solverSupported, "OperatorSmoother does not provide the smoother '" + type_ + "'.");
    if (solverSupported)
      SetParameterList(paramList);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
    Factory::SetParameterList(paramList);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");
    this->Input(currentLevel, "D0");
    this->Input(currentLevel, "M1");
    this->Input(currentLevel, "Ms");
    this->Input(currentLevel, "M0inv");
    this->Input(currentLevel, "Coordinates");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    SetupRefMaxwell(currentLevel);
    SmootherPrototype::IsSetup(true);
    this->GetOStream(Statistics1) << description() << std::endl;
  }


  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupRefMaxwell(Level& currentLevel) {

    using coordinateType = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
    using RealValuedMultiVector = Xpetra::MultiVector<coordinateType,LO,GO,NO>;

    auto SM        = currentLevel.Get<RCP<Matrix> >("A");
    auto D0        = currentLevel.Get<RCP<Matrix> >("D0");
    auto M1        = currentLevel.Get<RCP<Matrix> >("M1");
    auto Ms        = currentLevel.Get<RCP<Matrix> >("Ms");
    auto M0inv     = currentLevel.Get<RCP<Matrix> >("M0inv");
    auto coords    = currentLevel.Get<RCP<RealValuedMultiVector> >("Coordinates");
    auto params    = this->GetParameterList();

    op_ = rcp(new MueLu::RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>(SM, D0, Ms, M0inv, M1,
                                                                               Teuchos::null, coords,
                                                                               params));
    std::ostringstream oss;
    op_->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(oss)));
    cachedDescription_ = oss.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::OperatorSmoother::Apply(): Setup() has not been called");

    if (InitialGuessIsZero) {
      X.putScalar(0.0);
    }
    op_->apply(B, X);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    RCP<OperatorSmoother> smoother = rcp(new OperatorSmoother(*this) );
    smoother->SetParameterList(this->GetParameterList());
    return smoother;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    if (SmootherPrototype::IsSetup()) {
      out << op_->description();
      out << std::endl << std::endl;
      // out << cachedDescription_;
    } else {
      out << "Operator";
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
    }

    if (verbLevel & External) {
      if (op_ != Teuchos::null) {
        Teuchos::OSTab tab2(out);
        out << *op_ << std::endl << std::endl;
      }
    }

    // if (verbLevel & Debug) {
    //   out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
    //        << "-" << std::endl
    //        << "RCP<solver_>: " << tSolver_ << std::endl;
    // }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t OperatorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
    return Teuchos::OrdinalTraits<size_t>::invalid();
  }


} // namespace MueLu

#endif // MUELU_OPERATORSMOOTHER_DEF_HPP
