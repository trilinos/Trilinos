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
#ifndef MUELU_DIRECTSOLVER_DEF_HPP
#define MUELU_DIRECTSOLVER_DEF_HPP

#include <Xpetra_Utils.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_DirectSolver_decl.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"

#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_AmesosSmoother.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DirectSolver(const std::string& type, const Teuchos::ParameterList& paramListIn)
    : type_(type) {
    // The original idea behind all smoothers was to use prototype pattern. However, it does not fully work of the dependencies
    // calculation. Particularly, we need to propagate DeclareInput to proper prototypes. Therefore, both TrilinosSmoother and
    // DirectSolver do not follow the pattern exactly.
    // The difference is that in order to propagate the calculation of dependencies, we need to call a DeclareInput on a
    // constructed object (or objects, as we have two different code branches for Epetra and Tpetra). The only place where we
    // could construct these objects is the constructor. Thus, we need to store RCPs, and both TrilinosSmoother and DirectSolver
    // obtain a state: they contain RCP to smoother prototypes.
    sEpetra_ = Teuchos::null;
    sTpetra_ = Teuchos::null;

    ParameterList paramList = paramListIn;

    // We want DirectSolver to be able to work with both Epetra and Tpetra objects, therefore we try to construct both
    // Amesos and Amesos2 solver prototypes. The construction really depends on configuration options.
    bool triedEpetra = false, triedTpetra = false;
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
    sTpetra_ = rcp(new Amesos2Smoother(type_, paramList));
    TEUCHOS_TEST_FOR_EXCEPTION(sTpetra_.is_null(), Exceptions::RuntimeError, "Unable to construct Amesos2 direct solver");
    triedTpetra = true;
#endif
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)
    try {
      // GetAmesosSmoother masks the template argument matching, and simply throws if template arguments are incompatible with Epetra
      sEpetra_ = GetAmesosSmoother<SC,LO,GO,NO,LMO>(type_, paramList);
      TEUCHOS_TEST_FOR_EXCEPTION(sEpetra_.is_null(), Exceptions::RuntimeError, "Unable to construct Amesos direct solver");
    } catch (Exceptions::RuntimeError) {
      // AmesosSmoother throws if Scalar != double, LocalOrdinal != int, GlobalOrdinal != int
      this->GetOStream(Debug) << "Skipping AmesosSmoother construction due to incorrect type" << std::endl;
    }
    triedEpetra = true;
#endif

    // Check if we were able to construct at least one solver. In many cases that's all we need, for instance if a user
    // simply wants to use Tpetra only stack, never enables Amesos, and always runs Tpetra objects.
    TEUCHOS_TEST_FOR_EXCEPTION(!triedEpetra && !triedTpetra, Exceptions::RuntimeError, "Unable to construct direct solver. Plase enable (TPETRA and AMESOS2) or (EPETRA and AMESOS)");

    this->SetParameterList(paramList);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) {
    // We need to propagate SetFactory to proper place
    if (!sEpetra_.is_null()) sEpetra_->SetFactory(varName, factory);
    if (!sTpetra_.is_null()) sTpetra_->SetFactory(varName, factory);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    // Decide whether we are running in Epetra or Tpetra mode
    //
    // Theoretically, we could make this decision in the constructor, and create only
    // one of the smoothers. But we want to be able to reuse, so one can imagine a scenario
    // where one first runs hierarchy with tpetra matrix, and then with epetra.
    bool useTpetra = (currentLevel.lib() == Xpetra::UseTpetra);
    s_ = (useTpetra ? sTpetra_ : sEpetra_);
    TEUCHOS_TEST_FOR_EXCEPTION(s_.is_null(), Exceptions::RuntimeError, "Direct solver for " << (useTpetra ? "Tpetra" : "Epetra") << " was not constructed");

    s_->DeclareInput(currentLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::DirectSolver::Setup(): Setup() has already been called";

    s_->Setup(currentLevel);

    SmootherPrototype::IsSetup(true);

    this->SetParameterList(s_->GetParameterList());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");

    s_->Apply(X, B, InitialGuessIsZero);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    RCP<DirectSolver> newSmoo =  rcp(new DirectSolver(*this));

    // We need to be quite careful with Copy
    // We still want DirectSolver to follow Prototype Pattern, so we need to hide the fact that we do have some state
    if (!sEpetra_.is_null())
      newSmoo->sEpetra_ = sEpetra_->Copy();
    if (!sTpetra_.is_null())
      newSmoo->sTpetra_ = sTpetra_->Copy();

    // Copy the default mode
    newSmoo->s_ = (s_.get() == sTpetra_.get() ? newSmoo->sTpetra_ : newSmoo->sEpetra_);
    newSmoo->SetParameterList(this->GetParameterList());

    return newSmoo;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    if (s_ != Teuchos::null) {
      out << s_->description();
    } else {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab3(out);
      out << this->GetParameterList();
    }

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }

} // namespace MueLu

#endif // MUELU_DIRECTSOLVER_DEF_HPP
