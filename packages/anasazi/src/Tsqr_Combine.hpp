#ifndef __TSQR_Combine_hpp
#define __TSQR_Combine_hpp

#include <Tsqr_Blas.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_MatView.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class Combine
  /// \brief Operations for combining cache blocks or R factors in TSQR
  ///
  /// This class encapsulates four computational primitives for TSQR
  /// (in which R, R_1, and R_2 each represent an n x n upper
  /// triangular matrix, A represents an m x n cache block, and C_1
  /// and C_2 represent cache blocks with some number of columsn p):
  ///
  /// \li Factor [R; A] (factor_inner)
  /// \li Factor [R_1; R_2] (factor_pair)
  /// \li Apply Q factor of [R; A] to [C_1; C_2] (apply_inner)
  /// \li Apply Q factor of [R_1; R_2] to [C_1; C_2] (apply_pair)
  ///
  /// TODO (mfh 22 June 2010) Implement for the complex case.  Call
  /// {Z,C}GERC, etc., work with ScalarTraits< Scalar >::conj(tau[k])
  /// in apply_*(), etc.  Somehow I've had bad luck twice with
  /// reimplementing the float and double real arithmetic cases in
  /// C++, so I've left them in Fortran for now.  This class is just a
  /// wrapper.
  ///
  template< class Ordinal, class Scalar, bool is_complex = ScalarTraits< Scalar >::is_complex >
  class Combine {
  private:
    typedef MatView< Ordinal, Scalar > mat_view;
    typedef ConstMatView< Ordinal, Scalar > const_mat_view;

    BLAS< Ordinal, Scalar > blas_;
    LAPACK< Ordinal, Scalar > lapack_;

  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;

    Combine () {}

    /// Whether or not the QR factorizations computed by methods of
    /// this class produce an R factor with all nonnegative diagonal
    /// entries.  
    static bool QR_produces_R_factor_with_nonnegative_diagonal()
    { 
      return LAPACK< Ordinal, Scalar >::QR_produces_R_factor_with_nonnegative_diagonal; 
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
    /// \param trans [in] "N" means apply Q, anything else (such as
    ///   "T") means apply Q^T.  ("H", for applying Q^H, is not
    ///   currently supported.)
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
    apply_inner (const char* const trans,
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
		 Scalar work[]) const;

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
		  Scalar work[]) const;

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
		 Scalar work[]) const;
    
    /// Apply Q factor (or Q^T) of the 2*ncols_Q by ncols_Q matrix
    /// [R_top; R_bot] (stored in R_bot and tau) to the 2*ncols_Q by
    /// ncols_C matrix [C_top; C_bot].  The two blocks C_top and C_bot
    /// may have different leading dimensions (ldc_top resp. ldc_bot).
    void
    apply_pair (const char* const trans, 
		const Ordinal ncols_C, 
		const Ordinal ncols_Q, 
		const Scalar R_bot[], 
		const Ordinal ldr_bot,
		const Scalar tau[], 
		Scalar C_top[], 
		const Ordinal ldc_top, 
		Scalar C_bot[], 
		const Ordinal ldc_bot, 
		Scalar work[]) const;
  };

  // "Forward declaration"
  template< class Ordinal, class Scalar >
  class Combine< Ordinal, Scalar, false > {
  private:
    typedef MatView< Ordinal, Scalar > mat_view;
    typedef ConstMatView< Ordinal, Scalar > const_mat_view;

    BLAS< Ordinal, Scalar > blas_;
    LAPACK< Ordinal, Scalar > lapack_;

  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;

    Combine () {}

    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      return LAPACK< Ordinal, Scalar >::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    apply_inner (const char* const trans,
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
		 Scalar work[]) const;

    void
    factor_inner (const Ordinal m,
		  const Ordinal n,
		  Scalar R[],
		  const Ordinal ldr,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const;

    void
    factor_pair (const Ordinal n,
		 Scalar R_top[],
		 const Ordinal ldr_top,
		 Scalar R_bot[],
		 const Ordinal ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const;
    
    void
    apply_pair (const char* const trans, 
		const Ordinal ncols_C, 
		const Ordinal ncols_Q, 
		const Scalar R_bot[], 
		const Ordinal ldr_bot,
		const Scalar tau[], 
		Scalar C_top[], 
		const Ordinal ldc_top, 
		Scalar C_bot[], 
		const Ordinal ldc_bot, 
		Scalar work[]) const;
  };

} // namespace TSQR

#endif // __TSQR_Combine_hpp
