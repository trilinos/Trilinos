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
#ifndef MUELU_CONSTRAINT_DECL_HPP
#define MUELU_CONSTRAINT_DECL_HPP

#include <Teuchos_SerialDenseMatrix.hpp>

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Utilities_fwd.hpp"

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

  template<class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::sparseOps>
  class Constraint : public BaseClass {
#undef MUELU_CONSTRAINT_SHORT
#include "MueLu_UseShortNames.hpp"
  public:

  /*!
    @class Constraint class.
    @brief Class which contains the constraint space details
    */

    //! @name Setup methods.
    //@{

    /*! Setup constraint.
        \param B -- Fine nullspace vectors
        \param Bc -- Coarse nullspace vectors
        \param Ppattern -- Nonzero sparsity pattern for the prolongator
     */
    void Setup(const MultiVector& B, const MultiVector& Bc, RCP<const CrsGraph> Ppattern);

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
    RCP<MultiVector> X_;                                        //!< Overlapped coarse nullspace
    RCP<const CrsGraph> Ppattern_;                              //!< Nonzero sparsity pattern
    ArrayRCP<Teuchos::SerialDenseMatrix<LO,SC> > XXtInv_;       //!< Array storing \f$(Q_i Q_i^H)^{-1}\f$
  };

} // namespace MueLu

#define MUELU_CONSTRAINT_SHORT
#endif // MUELU_CONSTRAINT_DECL_HPP
