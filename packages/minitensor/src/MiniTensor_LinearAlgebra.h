// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_LinearAlgebra_h)
#define MiniTensor_LinearAlgebra_h

#include "MiniTensor_Tensor.h"

namespace minitensor {

///
/// Tensor Frobenius norm
/// \return \f$ \sqrt{A:A} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm(Tensor<T, N> const & A);

///
/// Tensor 1-norm
/// \return \f$ \max_{j \in {0,1,2}}\Sigma_{i=0}^2 |A_{ij}| \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_1(Tensor<T, N> const & A);

///
/// Tensor infinity-norm
/// \return \f$ \max_{i \in {0,1,2}}\Sigma_{j=0}^2 |A_{ij}| \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_infinity(Tensor<T, N> const & A);

///
/// 2nd-order tensor inverse
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse(Tensor<T, N> const & A);

///
/// 2nd-order tensor inverse using analitical expression for 2 and 3 dimensions
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse_fast23(Tensor<T, N> const & A);

///
/// 2nd-order tensor inverse using full pivoting, very accurate
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
inverse_full_pivot(Tensor<T, N> const & A);

///
/// Subtensor
/// \param i index
/// \param j index
/// \return Subtensor with i-row and j-col deleted.
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
subtensor(Tensor<T, N> const & A, Index const i, Index const j);

///
/// Swap row. Echange rows i and j in place
/// \param i index
/// \param j index
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
swap_row(Tensor<T, N> & A, Index const i, Index const j);

///
/// Swap column. Echange columns i and j in place
/// \param i index
/// \param j index
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
swap_col(Tensor<T, N> & A, Index const i, Index const j);

///
/// Determinant
/// \return \f$ \det A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
det(Tensor<T, N> const & A);

///
/// Trace
/// \return \f$ A:I \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
trace(Tensor<T, N> const & A);

///
/// First invariant, trace
/// \return \f$ I_A = A:I \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I1(Tensor<T, N> const & A);

///
/// Second invariant
/// \return \f$ II_A = \frac{1}{2}((I_A)^2-I_{A^2}) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I2(Tensor<T, N> const & A);

///
/// Third invariant
/// \return \f$ III_A = \det A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
I3(Tensor<T, N> const & A);

///
/// Exponential map.
/// \return \f$ \exp A \f$
///
template <typename T, Index N> Tensor<T, N> exp(Tensor<T, N> const &A);

///
/// Exponential map by Taylor series, radius of convergence is infinity
/// \return \f$ \exp A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
exp_taylor(Tensor<T, N> const & A);

///
/// Exponential map by squaring and scaling and Padé approximants.
/// See algorithm 10.20 in Functions of Matrices, N.J. Higham, SIAM, 2008.
/// \return \f$ \exp A \f$
///
template <typename T, Index N> Tensor<T, N> exp_pade(Tensor<T, N> const &A);

///
/// Logarithmic map.
/// \return \f$ \log A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log(Tensor<T, N> const & A);

///
/// Logarithmic map by Taylor series, converges for \f$ |A-I| < 1 \f$
/// \return \f$ \log A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_taylor(Tensor<T, N> const & A);

///
/// Logarithmic map by Gregory series,
/// converges for \f$ \min_i \text{Re} \lambda_i(A) > 0 \f$
/// \return \f$ \log A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_gregory(Tensor<T, N> const & A);

///
/// Logarithmic map for symmetric tensor.
/// \return \f$ \log A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_sym(Tensor<T, N> const & A);

///
/// Logarithmic map for symmetric tensor using eigenvalue decomposition.
/// \return \f$ \log A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_eig_sym(Tensor<T, N> const & A);

///
/// Logarithmic map of a rotation
/// \param R with \f$ R \in SO(3) \f$
/// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation(Tensor<T, N> const & R);

///
/// Logarithmic map of a 180 degree rotation
/// \param R with \f$ R \in SO(3) \f$
/// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation_pi(Tensor<T, N> const & R);

///
/// Apply Givens-Jacobi rotation on the left in place.
/// \param c and s for a rotation G in form [c, s; -s, c]
/// \param i and k indices for rows and columns where rotation is applied.
/// \param A tensor to rotate
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
givens_left(T const & c, T const & s, Index i, Index k, Tensor<T, N> & A);

///
/// Apply Givens-Jacobi rotation on the right in place.
/// \param c and s for a rotation G in form [c, s; -s, c]
/// \param i and k indices for rows and columns where rotation is applied.
/// \param A tensor to rotate
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
givens_right(T const & c, T const & s, Index i, Index k, Tensor<T, N> & A);

///
/// Apply rank-one update on the left in place
/// \f$ A = (I - beta v v^T) A \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
rank_one_left(T const & beta, Vector<T, N> const & v, Tensor<T, N> & A);

///
/// Apply rank-one update on the right in place
/// \f$ A = A (I - beta v v^T) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
rank_one_right(T const & beta, Vector<T, N> const & v, Tensor<T, N> & A);

///
/// Exponential map of a skew-symmetric tensor
/// \param r \f$ r \in so(3) \f$
/// \return \f$ R = \exp R \f$ with \f$ R \in SO(3) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
exp_skew_symmetric(Tensor<T, N> const & r);

///
/// Off-diagonal norm. Useful for SVD and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ \sqrt(\sum_i \sum_{j, j\neq i} a_{ij}^2) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
norm_off_diagonal(Tensor<T, N> const & A);

///
/// Arg max abs. Useful for inverse and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ (p,q) = arg max_{i,j} |a_{ij}| \f$
///
template <typename T, Index N>
std::pair<Index, Index> arg_max_abs(Tensor<T, N> const &A);

///
/// Arg max off-diagonal. Useful for SVD and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ (p,q) = arg max_{i \neq j} |a_{ij}| \f$
///
template <typename T, Index N>
std::pair<Index, Index> arg_max_off_diagonal(Tensor<T, N> const &A);

///
/// Sort and index. Useful for ordering singular values
/// and eigenvalues and corresponding vectors in the
/// respective decompositions.
/// \param u vector to sort
/// \return v P sorted vector, permutation matrix such that v = P^T u
///
template <typename T, Index N>
std::pair<Vector<T, N>, Tensor<T, N>> sort_permutation(Vector<T, N> const &u);

///
/// Singular value decomposition (SVD)
/// \return \f$ A = USV^T\f$
///
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>> svd(Tensor<T, N> const &A);

///
/// Project to O(N) (Orthogonal Group) using a Newton-type algorithm.
/// See Higham's Functions of Matrices p210 [2008]
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ R = \arg min_Q \|A - Q\|\f$
/// This algorithm projects a given tensor in GL(N) to O(N).
/// The rotation/reflection obtained through this projection is
/// the orthogonal component of the real polar decomposition
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
polar_rotation(Tensor<T, N> const & A);

///
/// Left polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_left(Tensor<T, N> const &A);

///
/// Right polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ RU = A \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_right(Tensor<T, N> const &A);

///
/// Left polar decomposition computed with eigenvalue decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_left_eig(Tensor<T, N> const &A);

///
/// R^3 right polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ RU = F \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> polar_right_eig(Tensor<T, N> const &A);

///
/// Left polar decomposition with matrix logarithm for V
/// \param F tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
///
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV(Tensor<T, N> const &F);

template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV_eig(Tensor<T, N> const &F);

///
/// Left polar decomposition with matrix logarithm for V using eig_spd_cos
/// \param F tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
///
template <typename T, Index N>
std::tuple<Tensor<T, N>, Tensor<T, N>, Tensor<T, N>>
polar_left_logV_lame(Tensor<T, N> const &F);

///
/// Logarithmic map using BCH expansion (4 terms)
/// \param v tensor
/// \param r tensor
/// \return Baker-Campbell-Hausdorff series up to 4 terms
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
bch(Tensor<T, N> const & v, Tensor<T, N> const & r);

///
/// Symmetric Schur algorithm for R^2.
/// \param \f$ A = [f, g; g, h] \in S(2) \f$
/// \return \f$ c, s \rightarrow [c, -s; s, c]\f$ diagonalizes A$
///
template <typename T>
std::pair<T, T> schur_sym(const T f, const T g, const T h);

///
/// Givens rotation. [c, -s; s, c] [a; b] = [r; 0]
/// \return c and s
///
template <typename T> std::pair<T, T> givens(T const &a, T const &b);

///
/// Eigenvalue decomposition for symmetric 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_sym(Tensor<T, N> const &A);

///
/// Eigenvalue decomposition for SPD 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_spd(Tensor<T, N> const &A);

///
/// Eigenvalue decomposition for SPD 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
/// This algorithm comes from the journal article
/// Scherzinger and Dohrmann, CMAME 197 (2008) 4007-4015
///
template <typename T, Index N>
std::pair<Tensor<T, N>, Tensor<T, N>> eig_spd_cos(Tensor<T, N> const &A);

///
/// Cholesky decomposition, rank-1 update algorithm
/// (Matrix Computations 3rd ed., Golub & Van Loan, p145)
/// \param A assumed symmetric tensor
/// \return G Cholesky factor A = GG^T and completed (bool)
/// algorithm ran to completion
///
template <typename T, Index N>
std::pair<Tensor<T, N>, bool> cholesky(Tensor<T, N> const &A);

///
/// Preconditioner types
///
enum class PreconditionerType
{
  UNDEFINED = 0,
  IDENTITY = 1,
  DIAGONAL = 2,
  MAX_ABS_ROW = 3,
};

///
/// Compute a preconditioner for improving the conditioning of a
/// linear system.
///
template <typename T, Index N, typename RHS>
std::pair<Tensor<T, N>, RHS> precon(PreconditionerType const pt,
                                    Tensor<T, N> const &A, RHS const &B);

///
/// Solve linear system of equations.
/// This is meant for the solution of small linear systems of equations
/// typically found in constitutive updates.
/// Right now the implementation is very inefficient (but accurate)
/// as it just uses the inverse function. It is intended to be used in
/// conjunction with Kokkos to take advantage of thread parallelism.
/// \param A assumed non-singular tensor
/// \param b rhs of the system Ax=b
/// \return x solution(s) to the system Ax=b
///
template <typename T, Index N, typename RHS>
RHS solve(Tensor<T, N> const &A, RHS const &b,
          PreconditionerType const pt = PreconditionerType::IDENTITY);

template<typename T, Index N, typename RHS>
KOKKOS_INLINE_FUNCTION
RHS
solve_full_pivot(Tensor<T, N> const & A, RHS const & b);

///
/// Condition number: ratio of largest to smalest singular values.
///
template <typename T, Index N> T cond(Tensor<T, N> const &A);

///
/// Reciprocal condition number: ratio of smallest to largest singular values.
///
template <typename T, Index N> T inv_cond(Tensor<T, N> const &A);

} // namespace minitensor

#include "MiniTensor_LinearAlgebra.i.h"
#include "MiniTensor_LinearAlgebra.t.h"

#endif // MiniTensor_LinearAlgebra_h
