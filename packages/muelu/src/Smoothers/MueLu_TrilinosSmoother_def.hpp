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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TrilinosSmoother(std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap, RCP<FactoryBase> AFact)
    : type_(type), paramList_(paramList), overlap_(overlap), AFact_(AFact)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(overlap_ < 0, Exceptions::RuntimeError, "overlap_ < 0");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get()); // TODO: also call Ifpack or Ifpack2::DeclareInput?
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //FactoryMonitor m(*this, "Setup Smoother");
    if (SmootherPrototype::IsSetup() == true) VerboseObject::GetOStream(Warnings0, 0) << "Warning: MueLu::DirectSolver::Setup(): Setup() has already been called";
    TEUCHOS_TEST_FOR_EXCEPTION(s_ != Teuchos::null, Exceptions::RuntimeError, "IsSetup() == false but s_ != Teuchos::null. This does not make sense");

    Xpetra::UnderlyingLib lib = currentLevel.Get< RCP<Matrix> >("A", AFact_.get())->getRowMap()->lib();

    if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)
      s_ = rcp( new Ifpack2Smoother(type_, paramList_, overlap_) );
      s_->SetFactory("A", AFact_);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external library availables for preconditionning Tpetra matrices. Compile MueLu with Ifpack2.");
#endif
    } else if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
      s_ = GetIfpackSmoother<SC,LO,GO,NO,LMO>(TrilinosSmoother::Ifpack2ToIfpack1Type(type_), TrilinosSmoother::Ifpack2ToIfpack1Param(paramList_), overlap_, AFact_);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external library availables for preconditionning Epetra matrices. Compile MueLu with Ifpack.");
#endif
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "lib != UseTpetra && lib != UseEpetra");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "");

    s_->Setup(currentLevel);

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");
    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "IsSetup() == true but s_ == Teuchos::null. This does not make sense");

    s_->Apply(X, B, InitialGuessIsZero);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new TrilinosSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Ifpack2ToIfpack1Type(std::string const & type) {
    if (type == "RELAXATION") { return "point relaxation stand-alone"; }
    if (type == "CHEBYSHEV")  { return "Chebyshev"; }
    if (type == "ILUT")        { return "ILU"; } //TODO: ILU is not a valid Ifpack2 type. This is just a temporary work-around to use TrilinosSmoother with Epetra + ILU

    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Cannot convert Ifpack2 preconditioner name to Ifpack: unkown type: " + type);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ParameterList TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Ifpack2ToIfpack1Param(Teuchos::ParameterList const & ifpack2List) {
    Teuchos::ParameterList ifpack1List = ifpack2List;

    if (ifpack2List.isParameter("relaxation: type")) {
      std::string relaxationType = ifpack2List.get<std::string>("relaxation: type");
      if (relaxationType == "Symmetric Gauss-Seidel") {
        ifpack1List.remove("relaxation: type");
        ifpack1List.set("relaxation: type", "symmetric Gauss-Seidel");
      }
    }

    return ifpack1List;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    //TODO
    //     if (s_ != Teuchos::null) {
    //       // Teuchos::OSTab tab2(out);
    //       s_->print(out, verbLevel);
    //     }

    //     if (verbLevel & Debug) {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }

    if (verbLevel & Parameters1) {
      out0 << "PrecType: " << type_ << std::endl;
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      out0 << "Overlap: " << overlap_ << std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "Epetra PrecType: " << Ifpack2ToIfpack1Type(type_) << std::endl
           << "Epetra Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << Ifpack2ToIfpack1Param(paramList_); }
    }
  }

} // namespace MueLu

#endif // MUELU_TRILINOSSMOOTHER_DEF_HPP
