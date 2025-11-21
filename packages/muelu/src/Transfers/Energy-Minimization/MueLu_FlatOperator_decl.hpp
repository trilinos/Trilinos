// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_FLATOPERATOR_DECL_HPP
#define MUELU_FLATOPERATOR_DECL_HPP

#include "MueLu_Constraint.hpp"
#include "Teuchos_RCP.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Operator.hpp"

#include "MueLu_VerboseObject.hpp"

namespace MueLu {

//! @brief Interprets a matrix as an operator that acts on a vector of nonzeros via SpGEMM.
/*!
  This class takes a matrix A and a constraint object.
  The constraint specifies a sparsity pattern.
  The FlatOperator applies A via

  Y = flat(A*operator(X))

  Here, operator(X) takes a vector X and interpretes it as a matrix with given sparsity pattern.
  Then a SpGEMM is performed.
  The entries corresponding to the sparsity pattern are extracted and returned in the vector Y.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class FlatOperator
  : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    public VerboseObject {
 public:
  //@{

  //! The Map associated with the domain of this operator, which must be compatible with X.getMap().
  virtual const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getDomainMap() const {
    return map_;
  }

  //! The Map associated with the range of this operator, which must be compatible with Y.getMap().
  virtual const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getRangeMap() const {
    return map_;
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
    return false;
  }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! A simple one-line description of this object.
  std::string description() const {
    return mat_->description();
  }

  //! Print the object with the given verbosity level to a FancyOStream.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    mat_->describe(out, verbLevel);
  }

  //@}

  //! @name Xpetra specific
  //@{

  FlatOperator() = default;

  //! TpetraOperator constructor to wrap a Tpetra::Operator object
  FlatOperator(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mat,
               const Teuchos::RCP<MueLu::Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>> constraint);

  void CheckMaps() {
  }

  void AllocateTemporaryMatrix() const {
    if (tempMat_.is_null()) {
      auto pattern = constraint_->GetPattern();
      tempMat_     = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(pattern);
      tempMat_->fillComplete();
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
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mat_;
  RCP<MueLu::Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>> constraint_;
  RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map_;
  mutable RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tempMat_;

};  // FlatOperator class

}  // namespace MueLu

#define MUELU_FLATOPERATOR_SHORT
#endif
