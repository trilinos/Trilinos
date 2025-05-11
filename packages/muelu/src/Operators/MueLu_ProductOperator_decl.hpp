// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PRODUCTOPERATOR_DECL_HPP
#define MUELU_PRODUCTOPERATOR_DECL_HPP

#include "Teuchos_RCP.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Operator.hpp"

#include "Xpetra_Utils.hpp"

namespace MueLu {

//! @brief Takes a sequence of operators and applies their product.
/*!
  Wrap operators op0, op1, op2, ... and apply modes mode0, mode1, mode2, ...
  into a single operator with apply
  Y : =(op0^mode0 * op1^mode1 * op2^mode2 * ...) X.
  */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class ProductOperator : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  //@{

  //! The Map associated with the domain of this operator, which must be compatible with X.getMap().
  virtual const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getDomainMap() const {
    return (modes_[ops_.size() - 1] == Teuchos::NO_TRANS) ? ops_[ops_.size() - 1]->getDomainMap() : ops_[ops_.size() - 1]->getRangeMap();
  }

  //! The Map associated with the range of this operator, which must be compatible with Y.getMap().
  virtual const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getRangeMap() const {
    return (modes_[0] == Teuchos::NO_TRANS) ? ops_[0]->getRangeMap() : ops_[0]->getDomainMap();
  }

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
      vary according to the values of \c alpha and \c beta. Specifically
      - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
      - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
   */
  virtual void
  apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
        Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
        Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  /// \brief Whether this operator supports applying the transpose or conjugate transpose.
  virtual bool hasTransposeApply() const {
    return true;
  }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! A simple one-line description of this object.
  std::string description() const {
    std::string descr("");
    for (auto it = ops_.begin(); it != ops_.end(); ++it) {
      descr += (*it)->description();
    }
    return descr;
  }

  //! Print the object with the given verbosity level to a FancyOStream.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    for (auto it = ops_.begin(); it != ops_.end(); ++it) {
      (*it)->describe(out, verbLevel);
    }
  }

  //@}

  //! @name Xpetra specific
  //@{

  ProductOperator() = default;

  //! TpetraOperator constructor to wrap a Tpetra::Operator object
  ProductOperator(std::vector<Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> ops,
                  std::vector<Teuchos::ETransp> modes)
    : ops_(ops)
    , modes_(modes) {
    TEUCHOS_ASSERT(ops_.size() >= 1);
    TEUCHOS_ASSERT(modes_.size() == ops_.size());

    CheckMaps();
    AllocateTemporaryMultiVectors(/*NumVectors=*/1);
  }

  ProductOperator(std::vector<Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> ops)
    : ops_(ops) {
    TEUCHOS_ASSERT(ops_.size() >= 1);

    for (auto it = ops_.begin(); it != ops_.end(); ++it) {
      modes_.push_back(Teuchos::NO_TRANS);
    }

    CheckMaps();
    AllocateTemporaryMultiVectors(/*NumVectors=*/1);
  }

  void CheckMaps() {
    for (size_t i = 0; i < ops_.size() - 1; ++i) {
      auto mapLeftOp  = (modes_[i] == Teuchos::NO_TRANS) ? ops_[i]->getDomainMap() : ops_[i]->getRangeMap();
      auto mapRightOp = (modes_[i + 1] == Teuchos::NO_TRANS) ? ops_[i + 1]->getRangeMap() : ops_[i + 1]->getDomainMap();
      TEUCHOS_ASSERT(mapLeftOp->isSameAs(*mapRightOp));
    }
  }

  void AllocateTemporaryMultiVectors(size_t NumVectors) const {
    if ((tempVecs_.size() == 0) || (tempVecs_[0]->getNumVectors() != NumVectors)) {
      tempVecs_.resize(0);
      for (size_t i = 0; i < ops_.size() - 1; ++i) {
        auto map = (modes_[i] == Teuchos::NO_TRANS) ? ops_[i]->getDomainMap() : ops_[i]->getRangeMap();
        tempVecs_.push_back(Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map, NumVectors));
      }
    }
  }

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &R) const {
    const auto one  = Teuchos::ScalarTraits<Scalar>::one();
    const auto zero = Teuchos::ScalarTraits<Scalar>::zero();

    apply(X, R, Teuchos::NO_TRANS, one, zero);
    R.update(one, B, -one);
  }

  //@}

 private:
  std::vector<Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> ops_;
  std::vector<Teuchos::ETransp> modes_;
  mutable std::vector<Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> tempVecs_;

};  // ProductOperator class

}  // namespace MueLu

#define MUELU_PRODUCTOPERATOR_SHORT
#endif
