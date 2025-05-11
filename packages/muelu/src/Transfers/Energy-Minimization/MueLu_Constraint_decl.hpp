// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CONSTRAINT_DECL_HPP
#define MUELU_CONSTRAINT_DECL_HPP

#include "Teuchos_ScalarTraits.hpp"

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_ProductOperator_fwd.hpp"

#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosXpetraAdapter.hpp>

namespace MueLu {

//! @brief Constraint space information for the potential prolongator
/*!
  This class implements an idea of the constrained space.  In energy
  minimization, constrained space is used simultaneously with the iterative
  method to construct the final prolongator. The space has two different
  constraints.

  \section pattern_constraint Nonzero pattern constraint

  Nonzero pattern constraint means that the final prolongator must have the
  provided nonzero pattern. This is achieved on each step of the iterative
  method by restricting the graph of the temporary prolongator to the desired
  pattern. It is implemented in the Apply function.  \note We do not update
  the graph of the provided temporary prolongator as this is a very expensive
  procedure. Rather, we extract its values and replace some of the values of
  the matrix with the correct graph.

  \section space_constraint Coarse space approximation constraint

  Generally, the coarse space constraint is presented by some matrix (X or Q)
  (see, for instance, the article by Mandel, Brezina and Vanek '99. It is
  well known that this matrix can be permuted to have a block diagonal form,
  where each block corresponds to a row in the prolongator. Specifically, let
  P be the prolongator, and Q be the constraint matrix. Then the constraint
  is generally written as \f$Q P = B,\f$ where B is the fine nullspace
  multivector. Q is a block diagonal matrix, \f$Q = diag(Q_1, ..., Q_n)\f$, where n
  is the number of rows in P. Each block Q_i is of size NSDim x nnz_i, where
  NSDim is the number of fine nullspace vectors, and nnz_i is the number of
  nonzero elements in the i-th row of P.

  To constrain the potential prolongator (with correct sparsity pattern, i.e.
  after the application of the nonzero pattern constraint), one updates its
  values as

              \f[P = P - Q^H(QQ^H)^{-1}QP.\f]

  Because of the block diagonal form of Q, this can be done row-by-row.
  \note For efficiency reasons, we store \f[(QQ^H)^{-1}\f] in the XXtInv_
  array. These matrices are dense, but have small size (NSDim x NSDim).
  */

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class Constraint
  : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    public BaseClass {
#undef MUELU_CONSTRAINT_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  /*!
    @class Constraint class.
    @brief Class which contains the constraint space details
    */

  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  Constraint() = default;

  virtual MagnitudeType ResidualNorm(const RCP<const Matrix> P) const = 0;

  //! The Map associated with the domain of this operator, which must be compatible with X.getMap().
  virtual const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getDomainMap() const {
    return X_->getDomainMap();
  }

  //! The Map associated with the range of this operator, which must be compatible with Y.getMap().
  virtual const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getRangeMap() const {
    return X_->getDomainMap();
  }

  //! @name Apply methods.
  //@{

  //! Apply constraint.
  virtual void apply(const MultiVector& P,
                     MultiVector& Projected,
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
                     Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R) const {
    const auto one  = Teuchos::ScalarTraits<Scalar>::one();
    const auto zero = Teuchos::ScalarTraits<Scalar>::zero();

    apply(X, R, Teuchos::NO_TRANS, one, zero);
    R.update(one, B, -one);
  }

  //@}

  RCP<const CrsGraph> GetPattern() const {
    return Ppattern_;
  }

  void SetPattern(RCP<const CrsGraph>& Ppattern) {
    Ppattern_ = Ppattern;
  }

  void SetX(RCP<Matrix>& X) {
    X_ = X;

    // Allocate memory
    temp1_ = MultiVectorFactory::Build(X_->getRangeMap(), 1);
    temp2_ = MultiVectorFactory::Build(X_->getRangeMap(), 1);
    temp3_ = MultiVectorFactory::Build(X_->getDomainMap(), 1);
  }

  RCP<Matrix> GetConstraintMatrix() {
    return X_;
  }

  void AssignMatrixEntriesToVector(const Matrix& P, const RCP<const CrsGraph>& pattern, MultiVector& vecP) const;

  void AssignMatrixEntriesToVector(const Matrix& P, MultiVector& vecP) const;

  RCP<Matrix> GetMatrixWithEntriesFromVector(MultiVector& vecP) const;

  void LeastSquaresSolve(const MultiVector& B, MultiVector& C) const;

 protected:
  void PrepareLeastSquaresSolve(bool singular = false);

 private:
  //! The constraints matrix
  RCP<Matrix> X_;

  //! Nonzero sparsity pattern
  RCP<const CrsGraph> Ppattern_;

  using MV = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using OP = Belos::OperatorT<MV>;
  RCP<Belos::LinearProblem<Scalar, MV, OP>> problem_;
  RCP<Belos::SolverManager<Scalar, MV, OP>> solver_;

  //! Prepare least-squares solve using Belos
  void PrepareLeastSquaresSolveBelos(bool singular);

  //! Perform least-squares solve using Belos
  void LeastSquaresSolveBelos(const MultiVector& B, MultiVector& C) const;

  //! Inverse of X*X^T
  RCP<Matrix> invXXt_;

  //! Prepare direct solution of least-squares problem
  void PrepareLeastSquaresSolveDirect(bool singular);

  //! Direct solve of least-squares problem
  void LeastSquaresSolveDirect(const MultiVector& B, MultiVector& C) const;

  // temporary memory
  RCP<MultiVector> temp1_, temp2_, temp3_;
};

}  // namespace MueLu

#define MUELU_CONSTRAINT_SHORT
#endif
