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

#include <Teuchos_SerialDenseMatrix.hpp>

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

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
class Constraint : public BaseClass {
#undef MUELU_CONSTRAINT_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  /*!
    @class Constraint class.
    @brief Class which contains the constraint space details
    */

  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  //! @name Setup methods.
  //@{

  /*! Setup constraint.
      \param B -- Fine nullspace vectors
      \param Bc -- Coarse nullspace vectors
      \param Ppattern -- Nonzero sparsity pattern for the prolongator
   */
  void Setup(const RCP<MultiVector>& B, const RCP<MultiVector>& Bc, RCP<const CrsGraph> Ppattern);

  MagnitudeType ResidualNorm(const RCP<const Matrix> P) const;

  //@}

  //! @name Apply methods.
  //@{

  //! Apply constraint.
  void Apply(const Matrix& P, Matrix& Projected) const;

  //@}

  RCP<const CrsGraph> GetPattern() const {
    return Ppattern_;
  }

 private:
  RCP<MultiVector> B_;
  RCP<MultiVector> Bc_;
  RCP<MultiVector> X_;                                    //!< Overlapped coarse nullspace
  RCP<const CrsGraph> Ppattern_;                          //!< Nonzero sparsity pattern
  ArrayRCP<Teuchos::SerialDenseMatrix<LO, SC> > XXtInv_;  //!< Array storing \f$(Q_i Q_i^H)^{-1}\f$
};

}  // namespace MueLu

#define MUELU_CONSTRAINT_SHORT
#endif  // MUELU_CONSTRAINT_DECL_HPP
