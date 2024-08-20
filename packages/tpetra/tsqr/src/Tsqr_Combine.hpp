// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_Combine.hpp
/// \brief Interface to TSQR's six computational kernels.

#ifndef TSQR_COMBINE_HPP
#define TSQR_COMBINE_HPP

#include "Tsqr_ApplyType.hpp"
#include "Tsqr_MatView.hpp"

namespace TSQR {
  /// \class Combine
  /// \brief Interface to TSQR's six computational kernels
  /// \author Mark Hoemmen
  ///
  /// This class provides the six computational primitives required by
  /// TSQR.  The primitives are as follows, in which R, R_1, and R_2
  /// each represent an n x n upper triangular matrix, A represents an
  /// m x n cache block, and C_1 and C_2 represent cache blocks with
  /// some number of columns p:
  ///
  /// <ul>
  /// <li> Factor A (factor_first) </li>
  /// <li> Apply Q factor of A to C (apply_first) </li>
  /// <li> Factor [R; A] (factor_inner) </li>
  /// <li> Factor [R_1; R_2] (factor_pair) </li>
  /// <li> Apply Q factor of [R; A] to [C_1; C_2] (apply_inner) </li>
  /// <li> Apply Q factor of [R_1; R_2] to [C_1; C_2] (apply_pair) </li>
  /// </ul>
  ///
  /// \tparam Ordinal Type of indices into matrices.
  /// \tparam Scalar Type of entries of matrices.
  ///
  /// TSQR includes two implementations of the Combine interface:
  ///
  /// <ul>
  /// <li> CombineDefault, which uses LAPACK and copies in and out of
  ///   scratch space that it owns, and </li>
  /// <li> CombineNative, a C++ in-place (no scratch space) generic
  ///   implementation) </li>
  /// </ul>
  ///
  /// There used to be a third implementation, CombineFortran, but it
  /// relied on a Fortran 9x compiler and was thus not often tested,
  /// so we removed it.
  template<class Ordinal, class Scalar>
  class Combine {
  public:
    //! Type of matrix entries.
    using scalar_type = Scalar;
    //! Type of (intraprocess) matrix indices.
    using ordinal_type = Ordinal;

    virtual ~Combine () = default;

    /// \brief Whether or not the QR factorizations computed by
    ///   methods of this class produce an R factor with all
    ///   nonnegative diagonal entries.
    virtual bool
    QR_produces_R_factor_with_nonnegative_diagonal () const = 0;

    /// \brief Best work array size.
    ///
    /// \param num_rows_Q [in] Number of rows in each block of the
    ///   matrix to factor.  ("Block" means the part of the matrix
    ///   passed directly to factor_first or factor_inner.)
    ///
    /// \param num_cols_Q [in] Number of columns of the matrix to
    ///   factor (the input/output matrix of factor_first or
    ///   factor_inner).
    ///
    /// \param num_cols_C [in] Number of columns of the matrix output
    ///   of apply_first, apply_inner, or apply_pair (use the max of
    ///   all three).
    virtual ordinal_type
    work_size (const ordinal_type num_rows_Q,
               const ordinal_type num_cols_Q,
               const ordinal_type num_cols_C) const = 0;

    /// \brief Factor the first cache block.
    ///
    /// Compute the QR factorization of the nrows by ncols matrix A
    /// (with leading dimension lda).  Overwrite the upper triangle of
    /// A with the resulting R factor, and the lower trapezoid of A
    /// (along with the length ncols tau array) with the implicitly
    /// stored Q factor.
    ///
    /// \param A [in/out] On input: the nrows by ncols matrix (in
    ///   column-major order, with leading dimension lda) to factor.
    ///   On output: upper triangle contains the R factor, and lower
    ///   part contains the implicitly stored Q factor.
    /// \param tau [out] Array of length ncols; on output, the
    ///   scaling factors for the Householder reflectors
    /// \param work [out] Workspace array of length ncols
    virtual void
    factor_first (const MatView<ordinal_type, Scalar>& A,
                  Scalar tau[],
                  Scalar work[],
                  const ordinal_type lwork) = 0;

    /// \brief Apply the result of factor_first() to C.
    ///
    /// Apply the Q factor, as computed by factor_first() and stored
    /// implicitly in A and tau, to the matrix C.
    virtual void
    apply_first (const ApplyType& applyType,
                 const MatView<ordinal_type, const Scalar>& A,
                 const Scalar tau[],
                 const MatView<ordinal_type, Scalar>& C,
                 Scalar work[],
                 const ordinal_type lwork) = 0;

    /// \brief Factor [R; A] for square upper triangular R and cache block A.
    ///
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
    /// \param R [inout] "Top" upper triangular n by n block $R_{k-1}$.
    ///   Overwritten with the new R factor $R_k$ of $[R_{k-1}; A_k]$.
    /// \param A [inout] "Bottom" dense m by n block $A_k$.  Overwritten
    ///   with the Householder reflectors representing the Q factor of
    ///   $[R_{k-1}; A_k]$.
    /// \param tau [out] Scaling factors of the Householder reflectors.
    ///   Corresponds to the TAU output of LAPACK's _GEQRF.
    /// \param work [out] Workspace (length >= n; don't need lwork or
    ///   workspace query)
    virtual void
    factor_inner (const MatView<ordinal_type, Scalar>& R,
                  const MatView<ordinal_type, Scalar>& A,
                  Scalar tau[],
                  Scalar work[],
                  const ordinal_type lwork) = 0;

    /// Apply the result of factor_inner().
    ///
    /// Apply the Q factor stored in [R; A] to [C_top; C_bot], where
    ///
    /// <ul>
    /// <li> A is     m       by ncols_Q, </li>
    /// <li> R is     ncols_Q by ncols Q, </li>
    /// <li> C_top is ncols_Q by ncols_C, and </li>
    /// <li> C_bot is m       by ncols_C. </li>
    /// </ul>
    ///
    /// The C blocks are allowed, but not required, to have different
    /// strides ("leading dimensions," in BLAS and LAPACK terms).  R
    /// is upper triangular, so we do not need an explicit version of
    /// R here.  The Householder reflectors representing the Q factor
    /// are stored compactly in A (specifically, in all of A, not just
    /// the lower triangle) and tau.
    ///
    /// \param apply_type [in] NoTranspose means apply Q, Transpose
    ///   means apply Q^T, and ConjugateTranspose means apply Q^H.
    /// \param A [in] m by ncols_Q matrix, in which the Householder
    ///   reflectors representing the Q factor are stored
    /// \param tau [in] array of length ncols_Q, storing the scaling
    ///   factors for the Householder reflectors representing Q
    /// \param C_top [inout]  ncols_Q by ncols_C matrix
    /// \param C_bot [inout]  m by ncols_C matrix
    /// \param work [out]     workspace array of length ncols_C
    virtual void
    apply_inner (const ApplyType& apply_type,
                 const MatView<ordinal_type, const Scalar>& A,
                 const Scalar tau[],
                 const MatView<ordinal_type, Scalar>& C_top,
                 const MatView<ordinal_type, Scalar>& C_bot,
                 Scalar work[],
                 const ordinal_type lwork) = 0;

    /// \brief Factor the pair of square upper triangular matrices
    ///   [R_top; R_bot].
    ///
    /// Store the resulting R factor in R_top, and the resulting
    /// Householder reflectors implicitly in R_bot and tau.
    virtual void
    factor_pair (const MatView<ordinal_type, Scalar>& R_top,
                 const MatView<ordinal_type, Scalar>& R_bot,
                 Scalar tau[],
                 Scalar work[],
                 const ordinal_type lwork) = 0;

    /// \brief Apply the result of \c factor_pair().
    ///
    /// Apply Q factor (or Q^T or Q^H) of the 2*ncols_Q by ncols_Q
    /// matrix [R_top; R_bot] (stored in R_bot and tau) to the
    /// 2*ncols_Q by ncols_C matrix [C_top; C_bot].  The two blocks
    /// C_top and C_bot need not be stored contiguously in memory, and
    /// they may have different strides ("leading dimensions," in BLAS
    /// and LAPACK terms).
    ///
    /// \param apply_type [in] NoTranspose means apply Q, Transpose
    ///   means apply Q^T, and ConjugateTranspose means apply Q^H.
    virtual void
    apply_pair (const ApplyType& apply_type,
                const MatView<ordinal_type, const Scalar>& R_bot,
                const Scalar tau[],
                const MatView<ordinal_type, Scalar>& C_top,
                const MatView<ordinal_type, Scalar>& C_bot,
                Scalar work[],
                const ordinal_type lwork) = 0;
  };

} // namespace TSQR

#endif // TSQR_COMBINE_HPP
