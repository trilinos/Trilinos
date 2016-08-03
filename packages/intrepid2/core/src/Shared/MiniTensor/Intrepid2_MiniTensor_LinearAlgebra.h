// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid2_MiniTensor_LinearAlgebra_h)
#define Intrepid2_MiniTensor_LinearAlgebra_h

#include "Intrepid2_MiniTensor_Tensor.h"

namespace Intrepid2 {

///
/// Tensor Frobenius norm
/// \return \f$ \sqrt{A:A} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm(Tensor<T, N, ES> const & A);

///
/// Tensor 1-norm
/// \return \f$ \max_{j \in {0,1,2}}\Sigma_{i=0}^2 |A_{ij}| \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm_1(Tensor<T, N, ES> const & A);

///
/// Tensor infinity-norm
/// \return \f$ \max_{i \in {0,1,2}}\Sigma_{j=0}^2 |A_{ij}| \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm_infinity(Tensor<T, N, ES> const & A);

///
/// 2nd-order tensor inverse
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N,ES>
inverse(Tensor<T, N, ES> const & A);

///
/// 2nd-order tensor inverse using analitical expression for 2 and 3 dimensions
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
inverse_fast23(Tensor<T, N, ES> const & A);

///
/// 2nd-order tensor inverse using full pivoting, very accurate
/// \param A nonsingular tensor
/// \return \f$ A^{-1} \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
inverse_full_pivot(Tensor<T, N, ES> const & A);

///
/// Subtensor
/// \param i index
/// \param j index
/// \return Subtensor with i-row and j-col deleted.
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
subtensor(Tensor<T, N, ES> const & A, Index const i, Index const j);

///
/// Swap row. Echange rows i and j in place
/// \param i index
/// \param j index
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
swap_row(Tensor<T, N, ES> & A, Index const i, Index const j);

///
/// Swap column. Echange columns i and j in place
/// \param i index
/// \param j index
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
swap_col(Tensor<T, N, ES> & A, Index const i, Index const j);

///
/// Determinant
/// \return \f$ \det A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
det(Tensor<T, N, ES> const & A);

///
/// Trace
/// \return \f$ A:I \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
trace(Tensor<T, N, ES> const & A);

///
/// First invariant, trace
/// \return \f$ I_A = A:I \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
I1(Tensor<T, N, ES> const & A);

///
/// Second invariant
/// \return \f$ II_A = \frac{1}{2}((I_A)^2-I_{A^2}) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
I2(Tensor<T, N, ES> const & A);

///
/// Third invariant
/// \return \f$ III_A = \det A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
I3(Tensor<T, N, ES> const & A);

///
/// Exponential map.
/// \return \f$ \exp A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
exp(Tensor<T, N, ES> const & A);

///
/// Exponential map by Taylor series, radius of convergence is infinity
/// \return \f$ \exp A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
exp_taylor(Tensor<T, N, ES> const & A);

///
/// Exponential map by squaring and scaling and Padé approximants.
/// See algorithm 10.20 in Functions of Matrices, N.J. Higham, SIAM, 2008.
/// \return \f$ \exp A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
exp_pade(Tensor<T, N, ES> const & A);

///
/// Logarithmic map.
/// \return \f$ \log A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
log(Tensor<T, N, ES> const & A);

///
/// Logarithmic map by Taylor series, converges for \f$ |A-I| < 1 \f$
/// \return \f$ \log A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
log_taylor(Tensor<T, N, ES> const & A);

///
/// Logarithmic map by Gregory series,
/// converges for \f$ \min_i \text{Re} \lambda_i(A) > 0 \f$
/// \return \f$ \log A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
log_gregory(Tensor<T, N, ES> const & A);

///
/// Logarithmic map for symmetric tensor.
/// \return \f$ \log A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
log_sym(Tensor<T, N, ES> const & A);

///
/// Logarithmic map for symmetric tensor using eigenvalue decomposition.
/// \return \f$ \log A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
log_eig_sym(Tensor<T, N, ES> const & A);

///
/// Logarithmic map of a rotation
/// \param R with \f$ R \in SO(3) \f$
/// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
log_rotation(Tensor<T, N, ES> const & R);

///
/// Logarithmic map of a 180 degree rotation
/// \param R with \f$ R \in SO(3) \f$
/// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
log_rotation_pi(Tensor<T, N, ES> const & R);

///
/// Apply Givens-Jacobi rotation on the left in place.
/// \param c and s for a rotation G in form [c, s; -s, c]
/// \param i and k indices for rows and columns where rotation is applied.
/// \param A tensor to rotate
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
givens_left(T const & c, T const & s, Index i, Index k, Tensor<T, N, ES> & A);

///
/// Apply Givens-Jacobi rotation on the right in place.
/// \param c and s for a rotation G in form [c, s; -s, c]
/// \param i and k indices for rows and columns where rotation is applied.
/// \param A tensor to rotate
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
givens_right(T const & c, T const & s, Index i, Index k, Tensor<T, N, ES> & A);

///
/// Apply rank-one update on the left in place
/// \f$ A = (I - beta v v^T) A \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
rank_one_left(T const & beta, Vector<T, N, ES> const & v, Tensor<T, N, ES> & A);

///
/// Apply rank-one update on the right in place
/// \f$ A = A (I - beta v v^T) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
void
rank_one_right(T const & beta, Vector<T, N, ES> const & v, Tensor<T, N, ES> & A);

///
/// Exponential map of a skew-symmetric tensor
/// \param r \f$ r \in so(3) \f$
/// \return \f$ R = \exp R \f$ with \f$ R \in SO(3) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
exp_skew_symmetric(Tensor<T, N, ES> const & r);

///
/// Off-diagonal norm. Useful for SVD and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ \sqrt(\sum_i \sum_{j, j\neq i} a_{ij}^2) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
T
norm_off_diagonal(Tensor<T, N, ES> const & A);

///
/// Arg max abs. Useful for inverse and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ (p,q) = arg max_{i,j} |a_{ij}| \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Index, Index>
arg_max_abs(Tensor<T, N, ES> const & A);

///
/// Arg max off-diagonal. Useful for SVD and other algorithms
/// that rely on Jacobi-type procedures.
/// \return \f$ (p,q) = arg max_{i \neq j} |a_{ij}| \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Index, Index>
arg_max_off_diagonal(Tensor<T, N, ES> const & A);

///
/// Sort and index. Useful for ordering singular values
/// and eigenvalues and corresponding vectors in the
/// respective decompositions.
/// \param u vector to sort
/// \return v P sorted vector, permutation matrix such that v = P^T u
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Vector<T, N, ES>, Tensor<T, N, ES>>
sort_permutation(Vector<T, N, ES> const & u);

///
/// Singular value decomposition (SVD)
/// \return \f$ A = USV^T\f$
///
template<typename T, Index N,  typename ES>
boost::tuple<Tensor<T, N, ES>, Tensor<T, N,ES>, Tensor<T, N,ES>>
svd(Tensor<T, N, ES> const & A);

///
/// Project to O(N) (Orthogonal Group) using a Newton-type algorithm.
/// See Higham's Functions of Matrices p210 [2008]
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ R = \arg min_Q \|A - Q\|\f$
/// This algorithm projects a given tensor in GL(N) to O(N).
/// The rotation/reflection obtained through this projection is
/// the orthogonal component of the real polar decomposition
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
polar_rotation(Tensor<T, N, ES> const & A);

///
/// Left polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, Tensor<T, N, ES>>
polar_left(Tensor<T, N, ES> const & A);

///
/// Right polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ RU = A \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, Tensor<T, N, ES>>
polar_right(Tensor<T, N, ES> const & A);

///
/// Left polar decomposition computed with eigenvalue decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, Tensor<T, N, ES>>
polar_left_eig(Tensor<T, N, ES> const & A);

///
/// R^3 right polar decomposition
/// \param A tensor (often a deformation-gradient-like tensor)
/// \return \f$ RU = F \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, Tensor<T, N, ES>>
polar_right_eig(Tensor<T, N, ES> const & A);

///
/// Left polar decomposition with matrix logarithm for V
/// \param F tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
///
template<typename T, Index N,  typename ES>
boost::tuple<Tensor<T, N, ES>, Tensor<T, N, ES>, Tensor<T, N, ES>>
polar_left_logV(Tensor<T, N, ES> const & F);

template<typename T, Index N,  typename ES>
boost::tuple<Tensor<T, N, ES>, Tensor<T, N, ES>, Tensor<T, N, ES>>
polar_left_logV_eig(Tensor<T, N, ES> const & F);

///
/// Left polar decomposition with matrix logarithm for V using eig_spd_cos
/// \param F tensor (often a deformation-gradient-like tensor)
/// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
///
template<typename T, Index N,  typename ES>
boost::tuple<Tensor<T, N, ES>, Tensor<T, N, ES>, Tensor<T, N, ES>>
polar_left_logV_lame(Tensor<T, N, ES> const & F);

///
/// Logarithmic map using BCH expansion (4 terms)
/// \param v tensor
/// \param r tensor
/// \return Baker-Campbell-Hausdorff series up to 4 terms
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Tensor<T, N, ES>
bch(Tensor<T, N, ES> const & v, Tensor<T, N, ES> const & r);

///
/// Symmetric Schur algorithm for R^2.
/// \param \f$ A = [f, g; g, h] \in S(2) \f$
/// \return \f$ c, s \rightarrow [c, -s; s, c]\f$ diagonalizes A$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
std::pair<T, T>
schur_sym(const T f, const T g, const T h);

///
/// Givens rotation. [c, -s; s, c] [a; b] = [r; 0]
/// \return c and s
///
template<typename T>
KOKKOS_INLINE_FUNCTION
std::pair<T, T>
givens(T const & a, T const & b);

///
/// Eigenvalue decomposition for symmetric 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, Tensor<T, N, ES>>
eig_sym(Tensor<T, N, ES> const & A);

///
/// Eigenvalue decomposition for SPD 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, Tensor<T, N, ES>>
eig_spd(Tensor<T, N, ES> const & A);

///
/// Eigenvalue decomposition for SPD 2nd-order tensor
/// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
/// This algorithm comes from the journal article
/// Scherzinger and Dohrmann, CMAME 197 (2008) 4007-4015
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, Tensor<T, N, ES>>
eig_spd_cos(Tensor<T, N, ES> const & A);

///
/// Cholesky decomposition, rank-1 update algorithm
/// (Matrix Computations 3rd ed., Golub & Van Loan, p145)
/// \param A assumed symmetric tensor
/// \return G Cholesky factor A = GG^T and completed (bool)
/// algorithm ran to completion
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
std::pair<Tensor<T, N, ES>, bool >
cholesky(Tensor<T, N, ES> const & A);

///
/// Solve linear system of equations.
/// This is means for the solution of small linear systems of equations
/// typically found in constitutive updates.
/// Right now the implementation is very inefficient (but accurate)
/// as it just uses the inverse function. It is intended to be used in
/// conjunction with Kokkos to take advantage of thread parallelism.
/// \param A assumed non-singular tensor
/// \param b rhs of the system Ax=b
/// \return x solution to the system Ax=b
///
template<typename T, Index N,  typename ES>
KOKKOS_INLINE_FUNCTION
Vector<T, N, ES>
solve(Tensor<T, N, ES> const & A, Vector<T, N, ES> const & b);

template<typename T, Index N, Index P,  typename ES>
KOKKOS_INLINE_FUNCTION
Matrix<T, N, P, ES>
solve(Tensor<T, N, ES> const & A, Matrix<T, N, P, ES> const & B);

///
template<typename T, Index N, typename ES>
KOKKOS_INLINE_FUNCTION
T
cond(Tensor<T, N, ES> const & A);

///
/// Reciprocal condition number: ratio of smallest to largest singular values.
///
template<typename T, Index N, typename ES>
KOKKOS_INLINE_FUNCTION
T
inv_cond(Tensor<T, N, ES> const & A);

} // namespace Intrepid2

#include "Intrepid2_MiniTensor_LinearAlgebra.i.h"
#include "Intrepid2_MiniTensor_LinearAlgebra.t.h"

#endif // Intrepid2_MiniTensor_LinearAlgebra_h
