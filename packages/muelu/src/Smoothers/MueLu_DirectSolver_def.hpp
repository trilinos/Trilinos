// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "MueLu_BelosSmoother.hpp"
#include "MueLu_StratimikosSmoother.hpp"
#include "MueLu_RefMaxwellSmoother.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
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
  sBelos_  = Teuchos::null;

  ParameterList paramList = paramListIn;

  // We want DirectSolver to be able to work with both Epetra and Tpetra objects, therefore we try to construct both
  // Amesos and Amesos2 solver prototypes. The construction really depends on configuration options.
  triedEpetra_ = triedTpetra_ = false;
#if defined(HAVE_MUELU_AMESOS2)
  try {
    sTpetra_ = rcp(new Amesos2Smoother(type_, paramList));
    if (sTpetra_.is_null())
      errorTpetra_ = "Unable to construct Amesos2 direct solver";
    else if (!sTpetra_->constructionSuccessful()) {
      errorTpetra_ = sTpetra_->constructionErrorMsg();
      sTpetra_     = Teuchos::null;
    }
  } catch (Exceptions::RuntimeError& e) {
    errorTpetra_ = e.what();
  } catch (Exceptions::BadCast& e) {
    errorTpetra_ = e.what();
  } catch (Teuchos::Exceptions::InvalidParameterName& e) {
    errorTpetra_ = e.what();
  }
  triedTpetra_ = true;
#endif
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)
  try {
    // GetAmesosSmoother masks the template argument matching, and simply throws if template arguments are incompatible with Epetra
    sEpetra_ = GetAmesosSmoother<SC, LO, GO, NO>(type_, paramList);
    if (sEpetra_.is_null())
      errorEpetra_ = "Unable to construct Amesos direct solver";
    else if (!sEpetra_->constructionSuccessful()) {
      errorEpetra_ = sEpetra_->constructionErrorMsg();
      sEpetra_     = Teuchos::null;
    }
  } catch (Exceptions::RuntimeError& e) {
    // AmesosSmoother throws if Scalar != double, LocalOrdinal != int, GlobalOrdinal != int
    errorEpetra_ = e.what();
  }
  triedEpetra_ = true;
#endif
#if defined(HAVE_MUELU_BELOS)
  try {
    sBelos_ = rcp(new BelosSmoother(type_, paramList));
    if (sBelos_.is_null())
      errorBelos_ = "Unable to construct Belos solver";
    else if (!sBelos_->constructionSuccessful()) {
      errorBelos_ = sBelos_->constructionErrorMsg();
      sBelos_     = Teuchos::null;
    }
  } catch (Exceptions::RuntimeError& e) {
    errorBelos_ = e.what();
  } catch (Exceptions::BadCast& e) {
    errorBelos_ = e.what();
  }
  triedBelos_ = true;
#endif
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  try {
    sStratimikos_ = rcp(new StratimikosSmoother(type_, paramList));
    if (sStratimikos_.is_null())
      errorStratimikos_ = "Unable to construct Stratimikos smoother";
    else if (!sStratimikos_->constructionSuccessful()) {
      errorStratimikos_ = sStratimikos_->constructionErrorMsg();
      sStratimikos_     = Teuchos::null;
    }
  } catch (Exceptions::RuntimeError& e) {
    errorStratimikos_ = e.what();
  }
  triedStratimikos_ = true;
#endif
  try {
    sRefMaxwell_ = rcp(new RefMaxwellSmoother(type_, paramList));
    if (sRefMaxwell_.is_null())
      errorRefMaxwell_ = "Unable to construct RefMaxwell smoother";
    else if (!sRefMaxwell_->constructionSuccessful()) {
      errorRefMaxwell_ = sRefMaxwell_->constructionErrorMsg();
      sRefMaxwell_     = Teuchos::null;
    }
  } catch (Exceptions::RuntimeError& e) {
    errorRefMaxwell_ = e.what();
  }
  triedRefMaxwell_ = true;

  // Check if we were able to construct at least one solver. In many cases that's all we need, for instance if a user
  // simply wants to use Tpetra only stack, never enables Amesos, and always runs Tpetra objects.
  TEUCHOS_TEST_FOR_EXCEPTION(!triedEpetra_ && !triedTpetra_ && !triedBelos_ && !triedStratimikos_ && !triedRefMaxwell_, Exceptions::RuntimeError,
                             "Unable to construct any direct solver."
                             "Plase enable (TPETRA and AMESOS2) or (EPETRA and AMESOS) or (BELOS) or (STRATIMIKOS)");

  TEUCHOS_TEST_FOR_EXCEPTION(sEpetra_.is_null() && sTpetra_.is_null() && sBelos_.is_null() && sStratimikos_.is_null() && sRefMaxwell_.is_null(), Exceptions::RuntimeError,
                             "Could not enable any direct solver:\n"
                                 << (triedEpetra_ ? "Epetra mode was disabled due to an error:\n" : "")
                                 << (triedEpetra_ ? errorEpetra_ : "")
                                 << (triedTpetra_ ? "Tpetra mode was disabled due to an error:\n" : "")
                                 << (triedTpetra_ ? errorTpetra_ : "")
                                 << (triedBelos_ ? "Belos was disabled due to an error:\n" : "")
                                 << (triedBelos_ ? errorBelos_ : "")
                                 << (triedStratimikos_ ? "Stratimikos was disabled due to an error:\n" : "")
                                 << (triedStratimikos_ ? errorStratimikos_ : "")
                                 << (triedRefMaxwell_ ? "RefMaxwell was disabled due to an error:\n" : "")
                                 << (triedRefMaxwell_ ? errorRefMaxwell_ : ""));
  ;

  this->SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) {
  // We need to propagate SetFactory to proper place
  if (!sEpetra_.is_null()) sEpetra_->SetFactory(varName, factory);
  if (!sTpetra_.is_null()) sTpetra_->SetFactory(varName, factory);
  if (!sBelos_.is_null()) sBelos_->SetFactory(varName, factory);
  if (!sStratimikos_.is_null()) sStratimikos_->SetFactory(varName, factory);
  if (!sRefMaxwell_.is_null()) sRefMaxwell_->SetFactory(varName, factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  if (!sBelos_.is_null())
    s_ = sBelos_;
  else if (!sStratimikos_.is_null())
    s_ = sStratimikos_;
  else if (!sRefMaxwell_.is_null())
    s_ = sRefMaxwell_;
  else {
    // Decide whether we are running in Epetra or Tpetra mode
    //
    // Theoretically, we could make this decision in the constructor, and create only
    // one of the smoothers. But we want to be able to reuse, so one can imagine a scenario
    // where one first runs hierarchy with tpetra matrix, and then with epetra.
    bool useTpetra = (currentLevel.lib() == Xpetra::UseTpetra);
    s_             = (useTpetra ? sTpetra_ : sEpetra_);
    if (s_.is_null()) {
      if (useTpetra) {
#if not defined(HAVE_MUELU_AMESOS2)
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                                   "Error: running in Tpetra mode, but MueLu with Amesos2 was disabled during the configure stage.\n"
                                   "Please make sure that:\n"
                                   "  - Amesos2 is enabled (Trilinos_ENABLE_Amesos2=ON),\n"
                                   "  - Amesos2 is available for MueLu to use (MueLu_ENABLE_Amesos2=ON)\n");
#else
        if (triedTpetra_)
          this->GetOStream(Errors) << "Tpetra mode was disabled due to an error:\n"
                                   << errorTpetra_ << std::endl;
#endif
      }
      if (!useTpetra) {
#if not defined(HAVE_MUELU_AMESOS)
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                                   "Error: running in Epetra mode, but MueLu with Amesos was disabled during the configure stage.\n"
                                   "Please make sure that:\n"
                                   "  - Amesos is enabled (you can do that with Trilinos_ENABLE_Amesos=ON),\n"
                                   "  - Amesos is available for MueLu to use (MueLu_ENABLE_Amesos=ON)\n");
#else
        if (triedEpetra_)
          this->GetOStream(Errors) << "Epetra mode was disabled due to an error:\n"
                                   << errorEpetra_ << std::endl;
#endif
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                                 "Direct solver for " << (useTpetra ? "Tpetra" : "Epetra") << " was not constructed");
    }
  }

  s_->DeclareInput(currentLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::DirectSolver::Setup(): Setup() has already been called" << std::endl;

  int oldRank = s_->SetProcRankVerbose(this->GetProcRankVerbose());

  s_->Setup(currentLevel);

  s_->SetProcRankVerbose(oldRank);

  SmootherPrototype::IsSetup(true);

  this->SetParameterList(s_->GetParameterList());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");

  s_->Apply(X, B, InitialGuessIsZero);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  RCP<DirectSolver> newSmoo = rcp(new DirectSolver(*this));

  // We need to be quite careful with Copy
  // We still want DirectSolver to follow Prototype Pattern, so we need to hide the fact that we do have some state
  if (!sEpetra_.is_null())
    newSmoo->sEpetra_ = sEpetra_->Copy();
  if (!sTpetra_.is_null())
    newSmoo->sTpetra_ = sTpetra_->Copy();
  if (!sBelos_.is_null())
    newSmoo->sBelos_ = sBelos_->Copy();
  if (!sStratimikos_.is_null())
    newSmoo->sStratimikos_ = sStratimikos_->Copy();
  if (!sRefMaxwell_.is_null())
    newSmoo->sRefMaxwell_ = sRefMaxwell_->Copy();

  // Copy the default mode
  if (s_.get() == sBelos_.get())
    newSmoo->s_ = newSmoo->sBelos_;
  else if (s_.get() == sStratimikos_.get())
    newSmoo->s_ = newSmoo->sStratimikos_;
  else if (s_.get() == sRefMaxwell_.get())
    newSmoo->s_ = newSmoo->sRefMaxwell_;
  else if (s_.get() == sTpetra_.get())
    newSmoo->s_ = newSmoo->sTpetra_;
  else
    newSmoo->s_ = newSmoo->sEpetra_;
  newSmoo->SetParameterList(this->GetParameterList());

  return newSmoo;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
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

}  // namespace MueLu

#endif  // MUELU_DIRECTSOLVER_DEF_HPP
