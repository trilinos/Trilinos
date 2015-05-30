/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

/// \file Ifpack2_Experimental_RBILUK_decl.hpp
/// \brief Declaration of RBILUK interface

#ifndef IFPACK2_EXPERIMENTALCRSRBILUK_DECL_HPP
#define IFPACK2_EXPERIMENTALCRSRBILUK_DECL_HPP

#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

#include <Ifpack2_RILUK.hpp>

namespace Ifpack2 {

namespace Experimental {

template <class scalar_type, class impl_scalar_type>
struct BlockMatrixOperations
{

  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  void square_matrix_matrix_multiply(const impl_scalar_type * a, const impl_scalar_type * b, impl_scalar_type * c,
      const int nrows, const impl_scalar_type alpha = STM::one(), const impl_scalar_type beta = STM::zero() ) const
  {
    for (int i = 0; i < nrows*nrows; ++i)
      c[i] = beta*c[i];

    for (int i = 0; i < nrows; ++i)
    {

      const int ioffset = i*nrows;
      for (int k = 0; k < nrows; ++k)
      {
        const int koffset = k*nrows;
        const impl_scalar_type val = alpha*a[ioffset+k];
        for (int j = 0; j < nrows; ++j)
        {
          c[ioffset+j] += val*b[koffset+j];
        }
      }
    }

  }


};

/** \class RBILUK
\brief ILU(k) factorization of a given Tpetra::Experimental::BlockCrsMatrix.
\tparam MatrixType A specialization of Tpetra::RowMatrix.

This class implements a "relaxed" incomplete ILU (ILU) factorization with level k fill.  It is based upon the ILU algorithms
outlined in Yousef Saad's "Iterative Methods for Sparse Linear Systems", 2nd edition, Chapter 10.

\section Ifpack2_Experimental_RBILUK_Parameters Parameters

For a complete list of valid parameters, see the documentation of setParameters().

The computed factorization is a function of several parameters:
<ul>
<li>
The graph structure (sparsity pattern) of the matrix: All fill is
derived from the original matrix nonzero structure.  Level zero fill
is defined as the original matrix pattern (nonzero structure), even if
the matrix value at an entry is stored as a zero. (Thus it is possible
to add entries to the ILU factors by adding zero entries to the
original matrix.)
</li>

<li>
Level of fill: Starting with the original matrix pattern as level
fill of zero, the next level of fill is determined by analyzing the
graph of the previous level and determining nonzero fill that is a
result of combining entries that were from previous level only (not
the current level).  This rule limits fill to entries that are direct
decendents from the previous level graph.  Fill for level k is
determined by applying this rule recursively.  For sufficiently large
values of k, the fill would eventually be complete and an exact LU
factorization would be computed.
</li>

<li>
Fraction of relaxation: Ifpack2::RILUK computes the ILU factorization
row-by-row.  As entries at a given row are computed, some number of
them will be dropped because they do match the prescribed sparsity
pattern.  The relaxation factor determines how these dropped values
will be handled.  If the factor is zero, then these extra entries will
by dropped.  This is a classical ILU approach.  If the RelaxValue is
1, then the sum of the extra entries will be added to the diagonal.
This is a classical Modified ILU (MILU) approach.  If RelaxValue is
between 0 and 1, then the factor times the sum of extra entries will
be added to the diagonal.

For most situations, the relaxation factor should be set to zero.  For
certain kinds of problems, e.g., reservoir modeling, there is a
conservation principle involved such that any operator should obey a
zero row-sum property.  MILU was designed for these cases and you
should set the relaxation factor to 1.  For other situations, setting
RelaxValue to some nonzero value may improve the stability of
factorization, and can be used if the computed ILU factors are poorly
conditioned.
</li>

<li>
Diagonal perturbation: Prior to computing the factorization, it is
possible to modify the diagonal entries of the matrix for which the
factorization will be computing.  If the absolute and relative
perturbation values are zero and one, respectively, the factorization
will be compute for the original user matrix A.  Otherwise, the
factorization will computed for a matrix that differs from the
original user matrix in the diagonal values only.  Below we discuss
the details of diagonal perturbations.
</li>

</ul>

*/
template<class MatrixType>
class RBILUK : virtual public Ifpack2::RILUK< Tpetra::RowMatrix< typename MatrixType::scalar_type,
  typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
{
 public:

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  typedef typename MatrixType::impl_scalar_type impl_scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The Node type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Tpetra::RowMatrix specialization used by this class.
  typedef Tpetra::RowMatrix<scalar_type,
      local_ordinal_type,
      global_ordinal_type,
      node_type> row_matrix_type;

  //! Tpetra::CrsMatrix specialization used by this class for representing L and U.
  typedef Tpetra::CrsMatrix<scalar_type,
      local_ordinal_type,
      global_ordinal_type,
      node_type> crs_matrix_type;

  typedef Tpetra::Experimental::BlockCrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> block_crs_matrix_type;

  template <class NewMatrixType> friend class RBILUK;

  /// \brief Constructor that takes a Tpetra::Experimental::BlockCrsMatrix
  ///
  /// \param A_in [in] The input matrix.
  RBILUK (const Teuchos::RCP<const block_crs_matrix_type>& A_in);

 private:
  /// \brief Copy constructor: declared private but not defined, so
  ///   that calling it is syntactically forbidden.
  RBILUK (const RBILUK<MatrixType> & src);

 public:

  //! Destructor (declared virtual for memory safety).
  virtual ~RBILUK ();

  //! Initialize by computing the symbolic incomplete factorization.
  void initialize ();

  /// \brief Compute the (numeric) incomplete factorization.
  ///
  /// This function computes the RBILU(k) factors L and U using the current:
  /// - Ifpack2_IlukGraph specifying the structure of L and U.
  /// - Value for the RILU(k) relaxation parameter.
  /// - Value for the a priori diagonal threshold values.
  ///
  /// initialize() must be called first, before this method may be called.
  void compute ();

  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  // Declare that we intend to overload RILUK::setMatrix, not hide it.
  // This avoids build warnings that the method below "hides
  // overloaded virtual function" (e.g., Clang 3.5).
  //
  // NOTE: If the base class of this class changes, e.g., if its
  // template parameter changes, then be sure to change the code below
  // to refer to the proper base class.
  using RILUK<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                typename MatrixType::local_ordinal_type,
                                typename MatrixType::global_ordinal_type,
                                typename MatrixType::node_type> >::setMatrix;

  /// \brief Change the matrix to be preconditioned.
  ///
  /// \param A [in] The new matrix.
  ///
  /// \post <tt>! isInitialized ()</tt>
  /// \post <tt>! isComputed ()</tt>
  ///
  /// Calling this method resets the preconditioner's state.  After
  /// calling this method with a nonnull input, you must first call
  /// initialize() and compute() (in that order) before you may call
  /// apply().
  ///
  /// You may call this method with a null input.  If A is null, then
  /// you may not call initialize() or compute() until you first call
  /// this method again with a nonnull input.  This method invalidates
  /// any previous factorization whether or not A is null, so calling
  /// setMatrix() with a null input is one way to clear the
  /// preconditioner's state (and free any memory that it may be
  /// using).
  ///
  /// The new matrix A need not necessarily have the same Maps or even
  /// the same communicator as the original matrix.
  void
  setMatrix (const Teuchos::RCP<const block_crs_matrix_type>& A);

  //@}
  //! @name Implementation of Teuchos::Describable interface
  //@{

  //! A one-line description of this object.
  std::string description () const;

  //@}
  //! \name Implementation of Tpetra::Operator
  //@{

  /// \brief Apply the (inverse of the) incomplete factorization to X, resulting in Y.
  ///
  /// For an incomplete factorization \f$A \approx LDU\f$, this method
  /// computes the following, depending on the value of \c mode:
  /// <ul>
  /// <li> If mode = Teuchos::NO_TRANS, it computes
  ///      <tt>Y = beta*Y + alpha*(U \ (D \ (L \ X)))</tt> </li>
  /// <li> If mode = Teuchos::TRANS, it computes
  ///      <tt>Y = beta*Y + alpha*(L^T \ (D^T \ (U^T \ X)))</tt> </li>
  /// <li> If mode = Teuchos::CONJ_TRANS, it computes
  ///      <tt>Y = beta*Y + alpha*(L^* \ (D^* \ (U^* \ X)))</tt>,
  ///      where the asterisk indicates the conjugate transpose. </li>
  /// </ul>
  /// If alpha is zero, then the result of applying the operator to a
  /// vector is ignored.  This matters because zero times NaN (not a
  /// number) is NaN, not zero.  Analogously, if beta is zero, then
  /// any values in Y on input are ignored.
  ///
  /// \param X [in] The input multivector.
  ///
  /// \param Y [in/out] The output multivector.
  ///
  /// \param mode [in] If Teuchos::TRANS resp. Teuchos::CONJ_TRANS,
  ///   apply the transpose resp. conjugate transpose of the incomplete
  ///   factorization.  Otherwise, don't apply the tranpose.
  ///
  /// \param alpha [in] Scaling factor for the result of applying the preconditioner.
  ///
  /// \param beta [in] Scaling factor for the initial value of Y.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one (),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero ()) const;
  //@}

public:

  //! Get the input matrix.
  Teuchos::RCP<const block_crs_matrix_type> getBlockMatrix () const;

  //! Return the L factor of the ILU factorization.
  const block_crs_matrix_type& getLBlock () const;

  //! Return the diagonal entries of the ILU factorization.
  const block_crs_matrix_type& getDBlock () const;

  //! Return the U factor of the ILU factorization.
  const block_crs_matrix_type& getUBlock () const;

private:
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;
  typedef Teuchos::ScalarTraits<impl_scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef typename block_crs_matrix_type::little_block_type little_block_type;
  typedef typename block_crs_matrix_type::little_vec_type little_vec_type;

  void allocate_L_and_U_blocks();
  void initAllValues (const block_crs_matrix_type& A);

  //! The (original) input matrix for which to compute ILU(k).
  Teuchos::RCP<const block_crs_matrix_type> A_block_;

  /// \brief The matrix used to to compute ILU(k).
  ///
  /// If A_local (the local filter of the original input matrix) is a
  /// Tpetra::CrsMatrix, then this is just A_local.  Otherwise, this
  /// class reserves the right for A_local_crs_ to be a copy of
  /// A_local.  This is because the current implementation of ILU(k)
  /// only knows how to factor a Tpetra::CrsMatrix.  That may change
  /// in the future.
  Teuchos::RCP<const block_crs_matrix_type> A_local_block_crs_;


  //! The L (lower triangular) factor of ILU(k).
  Teuchos::RCP<block_crs_matrix_type> L_block_;
  //! The U (upper triangular) factor of ILU(k).
  Teuchos::RCP<block_crs_matrix_type> U_block_;
  //! The diagonal entries of the ILU(k) factorization.
  Teuchos::RCP<block_crs_matrix_type> D_block_;

  //! The inverse of the diagonal
  Teuchos::RCP<block_crs_matrix_type> D_block_inverse_;

  BlockMatrixOperations<scalar_type,impl_scalar_type> blockMatOpts;

  void square_matrix_matrix_multiply(const impl_scalar_type * a, const impl_scalar_type * b, impl_scalar_type * c,
      const int nrows, const impl_scalar_type alpha = STM::one(), const impl_scalar_type beta = STM::zero() ) const
  {
    for (int i = 0; i < nrows*nrows; ++i)
      c[i] = beta*c[i];

    for (int i = 0; i < nrows; ++i)
    {

      const int ioffset = i*nrows;
      for (int k = 0; k < nrows; ++k)
      {
        const int koffset = k*nrows;
        const impl_scalar_type val = alpha*a[ioffset+k];
        for (int j = 0; j < nrows; ++j)
        {
          c[ioffset+j] += val*b[koffset+j];
        }
      }
    }

  }

};


} // namepsace Experimental

} // namespace Ifpack2

#endif /* IFPACK2_EXPERIMENTALCRSRBILUK_DECL_HPP */
