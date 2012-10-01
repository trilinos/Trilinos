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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DirectSolver(std::string const & type, Teuchos::ParameterList const & paramList, RCP<FactoryBase> AFact)
    : type_(type), paramList_(paramList), AFact_(AFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get()); // TODO: also call Amesos or Amesos2::DeclareInput?
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //FactoryMonitor m(*this, "Setup Smoother");
    if (SmootherPrototype::IsSetup() == true) VerboseObject::GetOStream(Warnings0, 0) << "Warning: MueLu::DirectSolver::Setup(): Setup() has already been called";
    TEUCHOS_TEST_FOR_EXCEPTION(s_ != Teuchos::null, Exceptions::RuntimeError, "IsSetup() == false but s_ != Teuchos::null. This does not make sense");

    Xpetra::UnderlyingLib lib = currentLevel.Get< RCP<Matrix> >("A", AFact_.get())->getRowMap()->lib();

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_AMESOS2
      s_ = rcp( new Amesos2Smoother(type_, paramList_, AFact_) );
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Tpetra matrices. Compile MueLu with Amesos2");
#endif
    } else if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_AMESOS
      s_ = GetAmesosSmoother<SC,LO,GO,NO,LMO>(type_, paramList_, AFact_);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Epetra matrices. Compile MueLu with Amesos"); // add Amesos2 to the msg when support for Amesos2+Epetra is implemented.
#endif
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "lib != UseTpetra && lib != UseEpetra");
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "");

    s_->Setup(currentLevel);

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");
    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "IsSetup() == true but s_ == Teuchos::null. This does not make sense");

    s_->Apply(X, B, InitialGuessIsZero);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new DirectSolver(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }
    
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
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
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab3(out); out << paramList_; }
    }

    if (verbLevel & Debug) {    
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu

#endif // MUELU_DIRECTSOLVER_DEF_HPP
