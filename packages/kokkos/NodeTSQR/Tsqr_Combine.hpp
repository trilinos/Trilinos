//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_Combine_hpp
#define __TSQR_Combine_hpp

#include <Teuchos_ScalarTraits.hpp>
#include <Tsqr_ApplyType.hpp>
#include <Tsqr_CombineNative.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class Combine
  /// \brief TSQR's six computational kernels
  ///
  /// This class encapsulates six computational primitives for TSQR
  /// (in which R, R_1, and R_2 each represent an n x n upper
  /// triangular matrix, A represents an m x n cache block, and C_1
  /// and C_2 represent cache blocks with some number of columns p):
  ///
  /// \li Factor A (factor_first)
  /// \li Apply Q factor of A to C (apply_first)
  /// \li Factor [R; A] (factor_inner)
  /// \li Factor [R_1; R_2] (factor_pair)
  /// \li Apply Q factor of [R; A] to [C_1; C_2] (apply_inner)
  /// \li Apply Q factor of [R_1; R_2] to [C_1; C_2] (apply_pair)
  ///
  /// CombineImpl is a particular implementation.  Its interface
  /// includes the same six public methods as does Combine's
  /// interface.  Combine itself just uses CombineImpl's
  /// implementations of these six methods.
  ///
  template< class Ordinal, 
	    class Scalar, 
	    class CombineImpl = CombineNative< Ordinal, Scalar, Teuchos::ScalarTraits< Scalar >::isComplex > >
  class Combine {
  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;
    typedef CombineImpl combine_impl_type;

    Combine () {}

    /// Whether or not the QR factorizations computed by methods of
    /// this class produce an R factor with all nonnegative diagonal
    /// entries.  
    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      return combine_impl_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    /// Compute the QR factorization of the nrows by ncols matrix A
    /// (with leading dimension lda).  Overwrite the upper triangle of
    /// A with the resulting R factor, and the lower trapezoid of A
    /// (along with the length ncols tau array) with the implicitly
    /// stored Q factor.  
    ///
    /// \param nrows [in] Number of rows in A
    /// \param ncols [in] Number of columns in A
    /// \param A [in/out] On input: the nrows by ncols matrix (in
    ///   column-major order, with leading dimension lda) to factor.
    ///   On output: upper triangle contains the R factor, and lower
    ///   part contains the implicitly stored Q factor.
    /// \param lda [in] Leading dimension of A
    /// \param tau [out] Array of length ncols; on output, the 
    ///   scaling factors for the Householder reflectors
    /// \param work [out] Workspace array of length ncols
    void
    factor_first (const Ordinal nrows,
		  const Ordinal ncols,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return impl_.factor_first (nrows, ncols, A, lda, tau, work);
    }

    /// Apply the Q factor, as computed by factor_first() and stored
    /// implicitly in A and tau, to the matrix C.
    void
    apply_first (const ApplyType& applyType,
		 const Ordinal nrows,
		 const Ordinal ncols_C,
		 const Ordinal ncols_A,
		 const Scalar A[],
		 const Ordinal lda,
		 const Scalar tau[],
		 Scalar C[],
		 const Ordinal ldc,
		 Scalar work[]) const
    {
      return impl_.apply_first (applyType, nrows, ncols_C, ncols_A, 
				A, lda, tau, C, ldc, work);
    }

    /// Apply the Q factor stored in [R; A] to [C_top; C_bot].  The C
    /// blocks are allowed, but not required, to have different leading
    /// dimensions (ldc_top resp. ldc_bottom).  R is upper triangular, so
    /// we do not need it; the Householder reflectors representing the Q
    /// factor are stored compactly in A (specifically, in all of A, not
    /// just the lower triangle).
    ///
    /// In the "sequential under parallel" version of TSQR, this function
    /// belongs to the sequential part (i.e., operating on cache blocks on
    /// a single processor).
    ///
    /// \param apply_type [in] NoTranspose means apply Q, Transpose
    ///   means apply Q^T, and ConjugateTranspose means apply Q^H.
    /// \param m [in]         number of rows of A
    /// \param ncols_C [in]   number of columns of [C_top; C_bot]
    /// \param ncols_Q [in]   number of columns of [R; A]
    /// \param A [in] m by ncols_Q matrix, in which the Householder 
    ///   reflectors representing the Q factor are stored
    /// \param lda [in]       leading dimension of A
    /// \param tau [in] array of length ncols_Q, storing the scaling
    ///   factors for the Householder reflectors representing Q 
    /// \param C_top [inout]  ncols_Q by ncols_C matrix
    /// \param ldc_top [in]   leading dimension of C_top
    /// \param C_bot [inout]  m by ncols_C matrix
    /// \param ldc_bot [in]   leading dimension of C_bot
    /// \param work [out]     workspace array of length ncols_C
    void
    apply_inner (const ApplyType& apply_type,
		 const Ordinal m,
		 const Ordinal ncols_C,
		 const Ordinal ncols_Q,
		 const Scalar A[],
		 const Ordinal lda,
		 const Scalar tau[],
		 Scalar C_top[],
		 const Ordinal ldc_top,
		 Scalar C_bot[],
		 const Ordinal ldc_bot,
		 Scalar work[]) const
    {
      impl_.apply_inner (apply_type, m, ncols_C, ncols_Q, 
			 A, lda, tau, 
			 C_top, ldc_top, C_bot, ldc_bot, work);
    }

    /// Perform one "inner" QR factorization step of sequential / parallel
    /// TSQR.  (In either case, only one processor calls this routine.)
    ///
    /// In the "sequential under parallel" version of TSQR, this function
    /// belongs to the sequential part (i.e., operating on cache blocks on
    /// a single processor).  Only the first cache block $A_0$ is factored
    /// as $Q_0 R_0 = A_0$ (see tsqr_factor_first); subsequent cache blocks
    /// $A_k$ are factored using this routine, which combines them with 
    /// $R_{k-1}$.
    ///
    /// Here is the matrix to factor:
    /// \[
    /// \begin{pmatrix}
    /// R_{k-1} \\      % $A_{k-1}$ is $m_{k-1} \times n$ with $m_{k-1} \geq n$
    /// A_k     \\      % $m_k \times n$ with $m_k \geq n$
    /// \end{pmatrix}
    /// \]
    ///
    /// Since $R_{k-1}$ is n by n upper triangular, we can store the
    /// Householder reflectors (representing the Q factor of $[R_{k-1};
    /// A_k]$) entirely in $A_k$ (specifically, in all of $A_k$, not just
    /// below the diagonal).
    ///
    /// \param m [in] Number of rows in the "bottom" block to factor.  
    ///   The number of rows in the top block doesn't matter, given the 
    ///   assumptions above, as long as $m_{k-1} \geq n$.
    /// \param n [in] Number of columns (same in both blocks)
    /// \param R [inout] "Top" upper triangular n by n block $R_{k-1}$.
    ///   Overwritten with the new R factor $R_k$ of $[R_{k-1}; A_k]$.
    /// \param ldr [in] Leading dimension of R
    /// \param A [inout] "Bottom" dense m by n block $A_k$.  Overwritten 
    ///   with the Householder reflectors representing the Q factor of 
    ///   $[R_{k-1}; A_k]$.
    /// \param tau [out] Scaling factors of the Householder reflectors.  
    ///   Corresponds to the TAU output of LAPACK's _GEQRF.
    /// \param work [out] Workspace (length >= n; don't need lwork or 
    ///   workspace query)
    void
    factor_inner (const Ordinal m,
		  const Ordinal n,
		  Scalar R[],
		  const Ordinal ldr,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      impl_.factor_inner (m, n, R, ldr, A, lda, tau, work);
    }

    /// Compute QR factorization of [R_top; R_bot].  Store resulting R
    /// factor in R_top and Householder reflectors in R_bot.
    ///
    /// \param n [in]         Number of rows and columns of each of R_top and R_bot
    /// \param R_top [inout]  n by n upper triangular matrix 
    /// \param ldr_top [in]   Leading dimension of R_top
    /// \param R_bot [inout]  n by n upper triangular matrix 
    /// \param ldr_bot [in]   Leading dimension of R_bot
    /// \param tau [out]      Scaling factors for Householder reflectors
    /// \param work [out]     Workspace array (of length >= n)
    ///
    void
    factor_pair (const Ordinal n,
		 Scalar R_top[],
		 const Ordinal ldr_top,
		 Scalar R_bot[],
		 const Ordinal ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const
    {
      impl_.factor_pair (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);      
    }
    
    /// Apply Q factor (or Q^T or Q^H) of the 2*ncols_Q by ncols_Q
    /// matrix [R_top; R_bot] (stored in R_bot and tau) to the
    /// 2*ncols_Q by ncols_C matrix [C_top; C_bot].  The two blocks
    /// C_top and C_bot may have different leading dimensions (ldc_top
    /// resp. ldc_bot).
    ///
    /// \param apply_type [in] NoTranspose means apply Q, Transpose
    ///   means apply Q^T, and ConjugateTranspose means apply Q^H.
    void
    apply_pair (const ApplyType& apply_type,
		const Ordinal ncols_C, 
		const Ordinal ncols_Q, 
		const Scalar R_bot[], 
		const Ordinal ldr_bot,
		const Scalar tau[], 
		Scalar C_top[], 
		const Ordinal ldc_top, 
		Scalar C_bot[], 
		const Ordinal ldc_bot, 
		Scalar work[]) const
    {
      impl_.apply_pair (apply_type, ncols_C, ncols_Q, 
			R_bot, ldr_bot, tau, 
			C_top, ldc_top, C_bot, ldc_bot, work);
    }

  private:
    combine_impl_type impl_;
  };

} // namespace TSQR

#endif // __TSQR_Combine_hpp
