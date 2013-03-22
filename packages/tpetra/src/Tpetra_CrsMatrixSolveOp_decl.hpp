// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIXSOLVEOP_HPP
#define TPETRA_CRSMATRIXSOLVEOP_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_CrsMatrix.hpp"

/*! \file Tpetra_CrsMatrixSolveOp_decl.hpp 

    Declaration of the class Tpetra::CrsMatrixSolveOp and related non-member constructors.
 */

namespace Tpetra {

  /// \class CrsMatrixSolveOp
  /// \brief Wrap a CrsMatrix instance's triangular solve in an Operator.
  ///
  /// \tparam Scalar Same as the first template parameter of Operator.
  ///   The type of the entries of the MultiVector input and output of
  ///   apply().  Not necessarily the same as the first template
  ///   parameter of the CrsMatrix used to create this object.
  /// \tparam MatScalar Same as the first template parameter of
  ///   CrsMatrix.  The type of the entries of the sparse matrix.  Not
  ///   necessarily the same as the type of the entries of the
  ///   MultiVector input and output of apply().
  /// \tparam LocalOrdinal Same as the second template parameter of
  ///   CrsMatrix and Operator.
  /// \tparam GlobalOrdinal Same as the third template parameter of
  ///   CrsMatrix and Operator.
  /// \tparam Node Same as the fourth template parameter of CrsMatrix
  ///   and Operator.
  /// \tparam LocalMatOps Same as the fifth template parameter of
  ///   CrsMatrix.
  ///
  /// This class' apply() method does a "local" triangular solve.
  /// "Local" is in quotes because apply() does the same communication
  /// (Import and Export) operations that CrsMatrix's apply() method
  /// would do for a sparse matrix-vector multiply, but the triangular
  /// solve is restricted to each process' part of the data.  Thus, it
  /// is not a triangular solve of a fully distributed triangular
  /// matrix.
  ///
  /// Here are some situations where this operation is useful:
  /// - Your sparse matrix A only lives in one MPI process, and you
  ///   have a factorization of it (either complete or incomplete).
  /// - Domain decomposition, where each MPI process owns one subdomain
  /// - Coarse-grid solves in algebraic multigrid
  /// - Mixed-precision operations, where the type <tt>MatScalar</tt>
  ///   of entries in the matrix differs from the type <tt>Scalar</tt>
  ///   of entries in the MultiVector input and output of apply().
  template <class Scalar, 
	    class MatScalar = Scalar, 
	    class LocalOrdinal = int, 
	    class GlobalOrdinal = LocalOrdinal, 
	    class Node = Kokkos::DefaultNode::DefaultNodeType, 
	    class LocalMatOps = typename Kokkos::DefaultKernels<MatScalar,LocalOrdinal,Node>::SparseOps>
  class CrsMatrixSolveOp : public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
  public:
    //! @name Constructor and destructor
    //@{ 

    //! Constructor; takes a CrsMatrix to use for local triangular solves.
    CrsMatrixSolveOp (const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A);

    //! Destructor
    virtual ~CrsMatrixSolveOp();

    //@}
    //! @name Implementation of Operator
    //@{ 

    /// \brief Compute \f$Y = \beta Y + \alpha B X\f$, where \f$B X\f$
    ///   represents the result of the local triangular solve.
    void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, 
	   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
	   Teuchos::ETransp mode = Teuchos::NO_TRANS, 
	   Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(), 
	   Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

    //! Whether apply() can solve with the (conjugate) transpose of the matrix.
    bool hasTransposeApply() const;

    /// The domain Map of this operator.
    /// This is the range map of the underlying CrsMatrix.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

    /// The range Map of this operator.
    /// This is the domain Map of the underlying CrsMatrix.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

    //@}
  protected:
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    //! The underlying CrsMatrix.
    const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > matrix_;

    //! Cached temporary destination of Import operation in apply().
    mutable Teuchos::RCP<MV> importMV_;
    //! Cached temporary source of Export operation in apply().
    mutable Teuchos::RCP<MV> exportMV_;

    //! Do the non-transpose solve.
    void 
    applyNonTranspose (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, 
		       MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y) const;
    //! Do the transpose solve.
    void 
    applyTranspose (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, 
		    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y) const;
  };

  /*! \brief Non-member function to create CrsMatrixSolveOp

      \relatesalso CrsMatrixSolveOp
   */
  template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP< CrsMatrixSolveOp<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >
  createCrsMatrixSolveOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A);

} // end of namespace Tpetra

#endif // TPETRA_CRSMATRIXSOLVEOP_HPP
