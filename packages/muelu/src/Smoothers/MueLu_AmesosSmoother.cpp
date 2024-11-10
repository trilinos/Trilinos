// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <algorithm>

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)

#include <Epetra_LinearProblem.h>

#include <Amesos_config.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Node>
AmesosSmoother<Node>::AmesosSmoother(const std::string& type, const Teuchos::ParameterList& paramList)
  : type_(type) {
  this->SetParameterList(paramList);

  if (!type_.empty()) {
    // Transform string to "Abcde" notation
    std::transform(type_.begin(), type_.end(), type_.begin(), ::tolower);
    std::transform(type_.begin(), ++type_.begin(), type_.begin(), ::toupper);
  }
  if (type_ == "Amesos_klu") type_ = "Klu";
  if (type_ == "Klu2") type_ = "Klu";
  if (type_ == "Amesos_umfpack") type_ = "Umfpack";
  if (type_ == "Superlu_dist") type_ = "Superludist";
  if (type_ == "Amesos_mumps") type_ = "Mumps";

  // Try to come up with something availble
  // Order corresponds to our preference
  // TODO: It would be great is Amesos provides directly this kind of logic for us
  std::string oldtype = type_;
  if (type_ == "" || Amesos().Query(type_) == false) {
#if defined(HAVE_AMESOS_SUPERLU)
    type_ = "Superlu";
#elif defined(HAVE_AMESOS_KLU)
    type_ = "Klu";
#elif defined(HAVE_AMESOS_SUPERLUDIST)
    type_ = "Superludist";
#elif defined(HAVE_AMESOS_UMFPACK)
    type_ = "Umfpack";
#else
    this->declareConstructionOutcome(true, "Amesos has been compiled without SuperLU_DIST, SuperLU, Umfpack or Klu. By default, MueLu tries" +
                                               "to use one of these libraries. Amesos must be compiled with one of these solvers,  " +
                                               "or a valid Amesos solver has to be specified explicitly.");
    return;
#endif
    if (oldtype != "")
      this->GetOStream(Warnings0) << "MueLu::AmesosSmoother: \"" << oldtype << "\" is not available. Using \"" << type_ << "\" instead" << std::endl;
    else
      this->GetOStream(Runtime1) << "MueLu::AmesosSmoother: using \"" << type_ << "\"" << std::endl;
  }
  this->declareConstructionOutcome(false, "");
}

template <class Node>
void AmesosSmoother<Node>::DeclareInput(Level& currentLevel) const {
  this->Input(currentLevel, "A");
}

template <class Node>
void AmesosSmoother<Node>::Setup(Level& currentLevel) {
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);

  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::AmesosSmoother::Setup(): Setup() has already been called" << std::endl;

  A_ = Factory::Get<RCP<Matrix> >(currentLevel, "A");

  RCP<Epetra_CrsMatrix> epA = Utilities::Op2NonConstEpetraCrs(A_);
  linearProblem_            = rcp(new Epetra_LinearProblem());
  linearProblem_->SetOperator(epA.get());

  Amesos factory;
  prec_ = rcp(factory.Create(type_, *linearProblem_));
  TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Solver '" + type_ + "' not supported by Amesos");

  // set Reindex flag, if A is distributed with non-contiguous maps
  // unfortunately there is no reindex for Amesos2, yet. So, this only works for Epetra based problems
  if (A_->getRowMap()->isDistributed() == true && A_->getRowMap()->isContiguous() == false)
    const_cast<ParameterList&>(this->GetParameterList()).set("Reindex", true);

  const ParameterList& paramList = this->GetParameterList();
  RCP<ParameterList> precList    = this->RemoveFactoriesFromList(paramList);

  prec_->SetParameters(*precList);

  const_cast<ParameterList&>(paramList).setParameters(*precList);

  int r = prec_->NumericFactorization();
  TEUCHOS_TEST_FOR_EXCEPTION(r != 0, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Amesos solver returns value of " + Teuchos::toString(r) + " during NumericFactorization()");

  SmootherPrototype::IsSetup(true);
}

template <class Node>
void AmesosSmoother<Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");

  Epetra_MultiVector& epX       = Utilities::MV2NonConstEpetraMV(X);
  Epetra_MultiVector const& epB = Utilities::MV2EpetraMV(B);
  // Epetra_LinearProblem takes the right-hand side as a non-const pointer.
  // I think this const_cast is safe because Amesos won't modify the rhs.
  Epetra_MultiVector& nonconstB = const_cast<Epetra_MultiVector&>(epB);

  linearProblem_->SetLHS(&epX);
  linearProblem_->SetRHS(&nonconstB);

  prec_->Solve();

  // Don't keep pointers to our vectors in the Epetra_LinearProblem.
  linearProblem_->SetLHS(0);
  linearProblem_->SetRHS(0);
}

template <class Node>
RCP<MueLu::SmootherPrototype<double, int, int, Node> > AmesosSmoother<Node>::Copy() const {
  return rcp(new AmesosSmoother<Node>(*this));
}

template <class Node>
std::string AmesosSmoother<Node>::description() const {
  std::ostringstream out;
  out << SmootherPrototype::description();
  out << "{type = " << type_ << "}";
  return out.str();
}

// using MueLu::Describable::describe; // overloading, not hiding
template <class Node>
void AmesosSmoother<Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0)
    out0 << "Prec. type: " << type_ << std::endl;

  if (verbLevel & Parameters1) {
    out0 << "Parameter list: " << std::endl;
    Teuchos::OSTab tab2(out);
    out << this->GetParameterList();
  }

  if (verbLevel & External)
    if (prec_ != Teuchos::null) {
      prec_->PrintStatus();
      prec_->PrintTiming();
    }

  if (verbLevel & Debug) {
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
         << "-" << std::endl
         << "RCP<A_>: " << A_ << std::endl
         << "RCP<linearProblem__>: " << linearProblem_ << std::endl
         << "RCP<prec_>: " << prec_ << std::endl;
  }
}

template <class Node>
size_t AmesosSmoother<Node>::getNodeSmootherComplexity() const {
  // FIXME: This is a placeholder
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

// The AmesosSmoother is only templated on the Node, since it is an Epetra only object
// Therefore we do not need the full ETI instantiations as we do for the other MueLu
// objects which are instantiated on all template parameters.
#if defined(HAVE_MUELU_EPETRA)
template class MueLu::AmesosSmoother<Xpetra::EpetraNode>;
#endif

#endif  // HAVE_MUELU_EPETRA && HAVE_MUELU_AMESOS
