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
#ifndef MUELU_AMESOS2SMOOTHER_DEF_HPP
#define MUELU_AMESOS2SMOOTHER_DEF_HPP

#include <algorithm>

#include "MueLu_ConfigDefs.hpp"
#if defined (HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Xpetra_Matrix.hpp>

#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include "MueLu_Amesos2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Amesos2Smoother(const std::string& type, const Teuchos::ParameterList& paramList)
    : type_(type) {
    this->SetParameterList(paramList);

    if (!type_.empty()) {
      // Transform string to "Abcde" notation
      std::transform(type_.begin(),   type_.end(),   type_.begin(), ::tolower);
      std::transform(type_.begin(), ++type_.begin(), type_.begin(), ::toupper);
    }
    if (type_ == "Superlu_dist")
      type_ = "Superludist";

    // Try to come up with something availble
    // Order corresponds to our preference
    // TODO: It would be great is Amesos2 provides directly this kind of logic for us
    if (type_ == "" || Amesos2::query(type_) == false) {
      std::string oldtype = type_;
#if defined(HAVE_AMESOS2_SUPERLU)
      type_ = "Superlu";
#elif defined(HAVE_AMESOS2_KLU2)
      type_ = "Klu";
#elif defined(HAVE_AMESOS2_SUPERLUDIST)
      type_ = "Superludist";
#elif defined(HAVE_AMESOS2_BASKER)
      type_ = "Basker";
#else
      throw Exceptions::RuntimeError("Amesos2 has been compiled without SuperLU_DIST, SuperLU, Klu, or Basker. By default, MueLu tries"
                                     "to use one of these libraries. Amesos2 must be compiled with one of these solvers, "
                                     "or a valid Amesos2 solver has to be specified explicitly.");
#endif
      if (oldtype != "")
        this->GetOStream(Warnings0) << "MueLu::Amesos2Smoother: \"" << oldtype << "\" is not available. Using \"" << type_ << "\" instead" << std::endl;
      else
        this->GetOStream(Runtime1) << "MueLu::Amesos2Smoother: using \"" << type_ << "\"" << std::endl;
    }

    // Check the validity of the solver type parameter
    TEUCHOS_TEST_FOR_EXCEPTION(Amesos2::query(type_) == false, Exceptions::RuntimeError, "The Amesos2 library reported that the solver '" << type_ << "' is not available. "
                               "Amesos2 has been compiled without the support of this solver, or the solver name is misspelled.");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Amesos2Smoother() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::Amesos2Smoother::Setup(): Setup() has already been called" << std::endl;

    RCP<Matrix> A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    RCP<Tpetra_CrsMatrix> tA = Utils::Op2NonConstTpetraCrs(A_);

    prec_ = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(type_, tA);
    TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "Amesos2::create returns Teuchos::null");

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::Amesos2Smoother::Apply(): Setup() has not been called");

    RCP<Tpetra_MultiVector> tX = Utils::MV2NonConstTpetraMV2(X);
    MultiVector & BNonC = const_cast<MultiVector&>(B);
    RCP<Tpetra_MultiVector> tB = Utils::MV2NonConstTpetraMV2(BNonC);
    prec_->setX(tX);
    prec_->setB(tB);

    prec_->solve();

    prec_->setX(Teuchos::null);
    prec_->setB(Teuchos::null);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Copy() const
  {
    return rcp (new Amesos2Smoother (*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;

    if (SmootherPrototype::IsSetup() == true) {
      out << prec_->description();

    } else {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
    }

    if ((verbLevel & External) && prec_ != Teuchos::null) {
      Teuchos::OSTab tab2(out);
      out << *prec_ << std::endl;
    }

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
  }

} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
#endif // MUELU_AMESOS2SMOOTHER_DEF_HPP
