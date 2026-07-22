// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_XPETRAOPERATOR_DEF_HPP
#define MUELU_XPETRAOPERATOR_DEF_HPP

#include "MueLu_XpetraOperator_decl.hpp"
#include "MueLu_Behavior.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    XpetraOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H)
  : Hierarchy_(H) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~XpetraOperator() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDomainMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;

  RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");
  return A->getDomainMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getRangeMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;

  RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");
  return A->getRangeMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
          Teuchos::ETransp mode,
          Scalar /* alpha */,
          Scalar /* beta */) const {
  TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS, std::logic_error, "MueLu::XpetraOperator does not support applying the adjoint operator");
  try {
    if (Behavior::debug()) {
      using Matrix  = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");

      // X is supposed to live in the range map of the operator (const rhs = B)
      TEUCHOS_TEST_FOR_EXCEPTION(A->getRangeMap()->isSameAs(*(X.getMap())) == false, std::logic_error,
                                 "MueLu::XpetraOperator::apply: map of X is incompatible with range map of A");
      TEUCHOS_TEST_FOR_EXCEPTION(A->getDomainMap()->isSameAs(*(Y.getMap())) == false, std::logic_error,
                                 "MueLu::XpetraOperator::apply: map of Y is incompatible with domain map of A");
    }
    Hierarchy_->Iterate(X, Y, 1, true);
  } catch (std::exception& e) {
    // FIXME add message and rethrow
    std::cerr << "Caught an exception in MueLu::XpetraOperator::apply():" << std::endl
              << e.what() << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    hasTransposeApply() const {
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
             Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R) const {
  using STS = Teuchos::ScalarTraits<Scalar>;
  R.update(STS::one(), B, STS::zero());
  this->apply(X, R, Teuchos::NO_TRANS, -STS::one(), STS::one());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetHierarchy() const {
  return Hierarchy_;
}

}  // namespace MueLu

#endif  // MUELU_XPETRAOPERATOR_DEF_HPP
