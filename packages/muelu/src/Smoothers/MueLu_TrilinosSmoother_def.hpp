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
#ifndef MUELU_TRILINOSSMOOTHER_DEF_HPP
#define MUELU_TRILINOSSMOOTHER_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_TrilinosSmoother_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TrilinosSmoother(const std::string& type, const Teuchos::ParameterList& paramListIn, const LO& overlap)
    : type_(type), overlap_(overlap)
  {
    // The original idea behind all smoothers was to use prototype pattern. However, it does not fully work of the dependencies
    // calculation. Particularly, we need to propagate DeclareInput to proper prototypes. Therefore, both TrilinosSmoother and
    // DirectSolver do not follow the pattern exactly.
    // The difference is that in order to propagate the calculation of dependencies, we need to call a DeclareInput on a
    // constructed object (or objects, as we have two different code branches for Epetra and Tpetra). The only place where we
    // could construct these objects is the constructor. Thus, we need to store RCPs, and both TrilinosSmoother and DirectSolver
    // obtain a state: they contain RCP to smoother prototypes.
    sEpetra_ = Teuchos::null;
    sTpetra_ = Teuchos::null;

    TEUCHOS_TEST_FOR_EXCEPTION(overlap_ < 0, Exceptions::RuntimeError, "Overlap parameter is negative (" << overlap << ")");

    ParameterList paramList = paramListIn;


    // We want TrilinosSmoother to be able to work with both Epetra and Tpetra objects, therefore we try to construct both
    // Ifpack and Ifpack2 smoother prototypes. The construction really depends on configuration options.
    bool triedEpetra = false, triedTpetra = false;
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)
    sTpetra_ = rcp(new Ifpack2Smoother(type_, paramList, overlap_));
    TEUCHOS_TEST_FOR_EXCEPTION(sTpetra_.is_null(), Exceptions::RuntimeError, "Unable to construct Ifpack2 smoother");
    triedTpetra = true;
#endif
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
    try {
      // GetIfpackSmoother masks the template argument matching, and simply throws if template arguments are incompatible with Epetra
      sEpetra_ = GetIfpackSmoother<SC,LO,GO,NO>(TrilinosSmoother::Ifpack2ToIfpack1Type(type_), TrilinosSmoother::Ifpack2ToIfpack1Param(paramList), overlap_);
      TEUCHOS_TEST_FOR_EXCEPTION(sEpetra_.is_null(), Exceptions::RuntimeError, "Unable to construct Ifpack smoother");
    } catch (Exceptions::RuntimeError) {
      // IfpackSmoother throws if Scalar != double, LocalOrdinal != int, GlobalOrdinal != int
      this->GetOStream(Debug) << "Skipping IfpackSmoother construction due to incorrect type" << std::endl;
    }
    triedEpetra = true;
#endif

    // Check if we were able to construct at least one smoother. In many cases that's all we need, for instance if a user
    // simply wants to use Tpetra only stack, never enables Ifpack, and always runs Tpetra objects.
    TEUCHOS_TEST_FOR_EXCEPTION(!triedEpetra && !triedTpetra,      Exceptions::RuntimeError, "Unable to construct Ifpack/Ifpack2 smoother."
                               "Plase enable (TPETRA and IFPACK2) or (EPETRA and IFPACK)");

    this->SetParameterList(paramList);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) {
    // We need to propagate SetFactory to proper place
    if (!sEpetra_.is_null()) sEpetra_->SetFactory(varName, factory);
    if (!sTpetra_.is_null()) sTpetra_->SetFactory(varName, factory);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    // Decide whether we are running in Epetra or Tpetra mode
    //
    // Theoretically, we could make this decision in the constructor, and create only
    // one of the smoothers. But we want to be able to reuse, so one can imagine a scenario
    // where one first runs hierarchy with tpetra matrix, and then with epetra.
    bool useTpetra = (currentLevel.lib() == Xpetra::UseTpetra);
    s_ = (useTpetra ? sTpetra_ : sEpetra_);
    TEUCHOS_TEST_FOR_EXCEPTION(s_.is_null(), Exceptions::RuntimeError, "Smoother for " << (useTpetra ? "Tpetra" : "Epetra") << " was not constructed");

    s_->DeclareInput(currentLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::TrilinosSmoother::Setup(): Setup() has already been called";

    s_->Setup(currentLevel);

    SmootherPrototype::IsSetup(true);

    this->SetParameterList(s_->GetParameterList());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::TrilinosSmoother::Apply(): Setup() has not been called");

    s_->Apply(X, B, InitialGuessIsZero);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    RCP<TrilinosSmoother> newSmoo = rcp(new TrilinosSmoother(type_, this->GetParameterList(), overlap_));

    // We need to be quite careful with Copy
    // We still want TrilinosSmoother to follow Prototype Pattern, so we need to hide the fact that we do have some state
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
  std::string TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Ifpack2ToIfpack1Type(const std::string& type) {
    if (type == "RELAXATION") { return "point relaxation stand-alone"; }
    if (type == "CHEBYSHEV")  { return "Chebyshev";                    }
    if (type == "ILUT")       { return "ILUT";                         }
    if (type == "RILUK")      { return "ILU";                          }
    if (type == "ILU")        { return "ILU";                          }
    if (type == "Amesos")     { return "Amesos";                       }

    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Cannot convert Ifpack2 preconditioner name to Ifpack: unknown type: " + type);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ParameterList TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Ifpack2ToIfpack1Param(const Teuchos::ParameterList& ifpack2List) {
    Teuchos::ParameterList ifpack1List = ifpack2List;

    if (ifpack2List.isParameter("relaxation: type") && ifpack2List.get<std::string>("relaxation: type") == "Symmetric Gauss-Seidel")
      ifpack1List.set("relaxation: type", "symmetric Gauss-Seidel");

    if (ifpack2List.isParameter("fact: iluk level-of-fill")) {
      ifpack1List.remove("fact: iluk level-of-fill");
      ifpack1List.set("fact: level-of-fill", ifpack2List.get<int>("fact: iluk level-of-fill"));
    }

    return ifpack1List;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
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
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Parameters1) {
      out0 << "PrecType: " << type_ << std::endl;
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
      out0 << "Overlap: " << overlap_ << std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "Epetra PrecType: " << Ifpack2ToIfpack1Type(type_) << std::endl
           << "Epetra Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << Ifpack2ToIfpack1Param(this->GetParameterList());;
    }
  }

} // namespace MueLu

#endif // MUELU_TRILINOSSMOOTHER_DEF_HPP
