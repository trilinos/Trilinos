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
#ifndef MUELU_IFPACK2SMOOTHER_DEF_HPP
#define MUELU_IFPACK2SMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include <Teuchos_ParameterList.hpp>

#include "Ifpack2_Factory.hpp"

#include "MueLu_Ifpack2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Ifpack2Smoother(std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap)
    : type_(type), paramList_(paramList), overlap_(overlap)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~Ifpack2Smoother() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameters(Teuchos::ParameterList const & paramList) {
    paramList_ = paramList;

    if (SmootherPrototype::IsSetup()) {
      // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
      // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...

      Teuchos::ParameterList nonConstParamList = paramList; // because Ifpack SetParameters() input argument is not const...
      prec_->setParameters(nonConstParamList);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ParameterList const & Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetParameters() { return paramList_; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);
    if (this->IsSetup() == true) this->GetOStream(Warnings0, 0) << "Warning: MueLu::Ifpack2Smoother::Setup(): Setup() has already been called";

    RCP<Matrix> A = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    Scalar negone = Teuchos::ScalarTraits<Scalar>::one() * Teuchos::as<Scalar>(-1.0);

    if (type_ == "CHEBYSHEV") {
      bool useCached=false;
      Scalar lambdaMax = negone;
      if ( !paramList_.isParameter("chebyshev: max eigenvalue") ) {
        lambdaMax = A->GetMaxEigenvalueEstimate();
        if (lambdaMax != negone) {
          useCached=true;
          this->GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue (cached with matrix)" << " = " << lambdaMax << std::endl;
        }
      }
      if (!useCached) {
        lambdaMax = paramList_.get("chebyshev: max eigenvalue",negone);
        if (lambdaMax != negone) {
          this->GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue (set in parameter list)" << " = " << lambdaMax << std::endl;
        } else {
          lambdaMax = Utils::PowerMethod(*A,true,10,1e-4);
          A->SetMaxEigenvalueEstimate(lambdaMax);
          this->GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue (calculated with power method)" << " = " << lambdaMax << std::endl;
        }
      }
      paramList_.set("chebyshev: max eigenvalue", lambdaMax);
    }

    RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> > tpA = Utils::Op2NonConstTpetraCrs(A);
    prec_ = Ifpack2::Factory::create(type_, tpA, overlap_);

    prec_->setParameters(paramList_);
    prec_->initialize();
    prec_->compute();

    SmootherPrototype::IsSetup(true);

    this->GetOStream(Statistics1, 0) << description() << std::endl;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Apply(): Setup() has not been called");

    // Forward the InitialGuessIsZero option to Ifpack2
    //  TODO: It might be nice to switch back the internal
    //        "zero starting solution" option of the ifpack2 object prec_ to his
    //        initial value at the end but there is no way right now to get
    //        the current value of the "zero starting solution" in ifpack2.
    //        It's not really an issue, as prec_  can only be used by this method.
    Teuchos::ParameterList  paramList = paramList_;
    if (type_ == "CHEBYSHEV") {
      paramList.set("chebyshev: zero starting solution", InitialGuessIsZero);
    }
    else if (type_ == "RELAXATION") {
      paramList.set("relaxation: zero starting solution", InitialGuessIsZero);
    }
    else if (type_ == "KRYLOV") {
      paramList.set("krylov: zero starting solution", InitialGuessIsZero);
    }
    else if (type_ == "SCHWARZ") {
      int overlap;
      Ifpack2::getParameter(paramList, "schwarz: overlap level", overlap);
      if (InitialGuessIsZero == false && overlap > 0) {
        if (this->IsPrint(Warnings0, 0)) {
          static int warning_only_once=0;
          if ((warning_only_once++) == 0)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::Ifpack2Smoother::Apply(): Additive Schwarz with overlap has no provision for a nonzero initial guess." << std::endl;
        }
      }
      else {
	paramList.set("schwarz: zero starting solution", InitialGuessIsZero);
      }
    }
    else if (type_ == "ILUT") {
      if (InitialGuessIsZero == false) {
        if (this->IsPrint(Warnings0, 0)) {
          static int warning_only_once=0;
          if ((warning_only_once++) == 0)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::Ifpack2Smoother::Apply(): ILUT has no provision for a nonzero initial guess." << std::endl;
          // TODO: ILUT using correction equation should be implemented in ifpack2 directly
          //       I think that an option named "zero starting solution"
          //       is also appropriate for ILUT
        }
      }
    } else {
      // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
      // we should remove the if/else/elseif and just test if this
      // option is supported by current ifpack2 preconditioner
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,"Ifpack2Smoother::Apply(): Ifpack2 preconditioner '"+type_+"' not supported");
    }
    prec_->setParameters(paramList);

    // Apply
    Tpetra::MultiVector<SC,LO,GO,NO> &tpX = Utils::MV2NonConstTpetraMV(X);
    Tpetra::MultiVector<SC,LO,GO,NO> const &tpB = Utils::MV2TpetraMV(B);
    prec_->apply(tpB,tpX);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp(new Ifpack2Smoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    if (SmootherPrototype::IsSetup()) {
      out << prec_->description();
    } else {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      out0 << "Overlap: "        << overlap_ << std::endl;
    }

    if (verbLevel & External) {
      if (prec_ != Teuchos::null) { Teuchos::OSTab tab2(out); out << *prec_ << std::endl; }
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_IFPACK2
#endif // MUELU_IFPACK2SMOOTHER_DEF_HPP
