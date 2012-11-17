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
#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_AMESOS

#include <Epetra_LinearProblem.h>

#include <Amesos_config.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  AmesosSmoother::AmesosSmoother(std::string const & type, Teuchos::ParameterList const & paramList, RCP<FactoryBase> AFact)
    : type_(type), paramList_(paramList), AFact_(AFact)
  {
    // set default solver type
    if(type_ == "") {
#if defined(HAVE_AMESOS_SUPERLU)
      type_ = "Superlu";     // 1. default smoother (if Superlu is available)
#elif defined(HAVE_AMESOS_KLU)
      type_ = "Klu";         // 2. default smoother (if KLU is available)
#elif defined(HAVE_AMESOS_SUPERLUDIST)
      type_ = "Superludist"; // 3. default smoother (if Superludist is available)
#elif defined(HAVE_AMESOS_UMFPACK)
      type_ = "Umfpack";     // 4. default smoother (if Umfpack is available)
#endif
    } // if(type_ == "")

    // check for valid direct solver type
    TEUCHOS_TEST_FOR_EXCEPTION(type_ != "Superlu" && type_ != "Superludist" && type_ != "Klu" && type_ != "Amesos_Klu" && type_ != "Umfpack" && type_ != "Amesos_Umfpack", Exceptions::RuntimeError, "MueLu::AmesosSmoother::AmesosSmoother(): Solver '" + type_ + "' not supported");
    if (type_ == "Superlu") {
#if not defined(HAVE_AMESOS_SUPERLU)
      TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::AmesosSmoother(): Amesos compiled without SuperLU. Cannot define a solver by default for this AmesosSmoother object");
#endif
    }
    if (type_ == "Superludist") {
#if not defined(HAVE_AMESOS_SUPERLUDIST)
      TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::AmesosSmoother(): Amesos compiled without SuperLU_DIST. Cannot define a solver by default for this AmesosSmoother object");
#endif
    }
    if (type_ == "Klu" || type_ == "Amesos_Klu") {
#if not defined(HAVE_AMESOS_KLU)
      TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::AmesosSmoother(): Amesos compiled without KLU. Cannot define a solver by default for this AmesosSmoother object");
#endif
    }
    if (type_ == "Umfpack" || type_ == "Amesos_Umfpack") {
#if not defined(HAVE_AMESOS_UMFPACK)
      TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::AmesosSmoother(): Amesos compiled without Umfpack. Cannot define a solver by default for this AmesosSmoother object");
#endif
    }
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == true, Exceptions::RuntimeError, "TO BE REMOVED");
  }

  AmesosSmoother::~AmesosSmoother() {}

  void AmesosSmoother::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
  }

  void AmesosSmoother::Setup(Level &currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true) GetOStream(Warnings0, 0) << "Warning: MueLu::AmesosSmoother::Setup(): Setup() has already been called" << std::endl;

    A_ = currentLevel.Get< RCP<Matrix> >("A", AFact_.get());

    RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
    linearProblem_ = rcp( new Epetra_LinearProblem() );
    linearProblem_->SetOperator(epA.get());

    Amesos factory;
    prec_ = rcp(factory.Create(type_, *linearProblem_));
    TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Solver '" + type_ + "' not supported by Amesos");

    // set Reindex flag, if A is distributed with non-contiguous maps
    // unfortunately there is no reindex for Amesos2, yet. So, this only works for Epetra based problems
    if(A_->getRowMap()->isDistributed() == true && A_->getRowMap()->isContiguous() == false)
      paramList_.set("Reindex", true);

    prec_->SetParameters(paramList_);

    int r = prec_->NumericFactorization();
    TEUCHOS_TEST_FOR_EXCEPTION(r != 0, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Amesos solver returns value of " + Teuchos::Utils::toString(r) + " during NumericFactorization()");

    SmootherPrototype::IsSetup(true);
  }

  void AmesosSmoother::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");

    Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
    Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);
    //Epetra_LinearProblem takes the right-hand side as a non-const pointer.
    //I think this const_cast is safe because Amesos won't modify the rhs.
    Epetra_MultiVector &nonconstB = const_cast<Epetra_MultiVector&>(epB);

    linearProblem_->SetLHS(&epX);
    linearProblem_->SetRHS(&nonconstB);

    prec_->Solve();

    // Don't keep pointers to our vectors in the Epetra_LinearProblem.
    linearProblem_->SetLHS(0);
    linearProblem_->SetRHS(0);
  }

  RCP<MueLu::SmootherPrototype<double,int,int> > AmesosSmoother::Copy() const {
    return rcp( new AmesosSmoother(*this) );
  }

  std::string AmesosSmoother::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  //using MueLu::Describable::describe; // overloading, not hiding
  void AmesosSmoother::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
    }

    if (verbLevel & External) {
      if (prec_ != Teuchos::null) { prec_->PrintStatus(); prec_->PrintTiming(); } //TODO: redirect output?
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<A_>: " << A_ << std::endl
           << "RCP<linearProblem__>: " << linearProblem_ << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

} // namespace MueLu

#endif
