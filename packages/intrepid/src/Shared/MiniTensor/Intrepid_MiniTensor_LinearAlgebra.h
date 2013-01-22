// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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

#if !defined(Intrepid_MiniTensor_LinearAlgebra_h)
#define Intrepid_MiniTensor_LinearAlgebra_h

#include "Intrepid_MiniTensor_Tensor.h"

namespace Intrepid {

  ///
  /// Tensor Frobenius norm
  /// \return \f$ \sqrt{A:A} \f$
  ///
  template<typename T>
  T
  norm(Tensor<T> const & A);

  ///
  /// Tensor 1-norm
  /// \return \f$ \max_{j \in {0,1,2}}\Sigma_{i=0}^2 |A_{ij}| \f$
  ///
  template<typename T>
  T
  norm_1(Tensor<T> const & A);

  ///
  /// Tensor infinity-norm
  /// \return \f$ \max_{i \in {0,1,2}}\Sigma_{j=0}^2 |A_{ij}| \f$
  ///
  template<typename T>
  T
  norm_infinity(Tensor<T> const & A);

  ///
  /// 2nd-order tensor inverse
  /// \param A nonsingular tensor
  /// \return \f$ A^{-1} \f$
  ///
  template<typename T>
  Tensor<T>
  inverse(Tensor<T> const & A);

  ///
  /// Subtensor
  /// \param A tensor
  /// \param i index
  /// \param j index
  /// \return Subtensor with i-row and j-col deleted.
  ///
  template<typename T>
  Tensor<T>
  subtensor(Tensor<T> const & A, Index const i, Index const j);

  ///
  /// Swap row. Echange rows i and j in place
  /// \param A tensor
  /// \param i index
  /// \param j index
  ///
  template<typename T>
  void
  swap_row(Tensor<T> & A, Index const i, Index const j);

  ///
  /// Swap column. Echange columns i and j in place
  /// \param A tensor
  /// \param i index
  /// \param j index
  ///
  template<typename T>
  void
  swap_col(Tensor<T> & A, Index const i, Index const j);

  ///
  /// Determinant
  /// \param A tensor
  /// \return \f$ \det A \f$
  ///
  template<typename T>
  T
  det(Tensor<T> const & A);

  ///
  /// Trace
  /// \param A tensor
  /// \return \f$ A:I \f$
  ///
  template<typename T>
  T
  trace(Tensor<T> const & A);

  ///
  /// First invariant, trace
  /// \param A tensor
  /// \return \f$ I_A = A:I \f$
  ///
  template<typename T>
  T
  I1(Tensor<T> const & A);

  ///
  /// Second invariant
  /// \param A tensor
  /// \return \f$ II_A = \frac{1}{2}((I_A)^2-I_{A^2}) \f$
  ///
  template<typename T>
  T
  I2(Tensor<T> const & A);

  ///
  /// Third invariant
  /// \param A tensor
  /// \return \f$ III_A = \det A \f$
  ///
  template<typename T>
  T
  I3(Tensor<T> const & A);

  ///
  /// Exponential map.
  /// \param A tensor
  /// \return \f$ \exp A \f$
  ///
  template<typename T>
  Tensor<T>
  exp(Tensor<T> const & A);

  ///
  /// Exponential map by Taylor series, radius of convergence is infinity
  /// \param A tensor
  /// \return \f$ \exp A \f$
  ///
  template<typename T>
  Tensor<T>
  exp_taylor(Tensor<T> const & A);

  ///
  /// Exponential map by squaring and scaling and Padé approximants.
  /// See algorithm 10.20 in Functions of Matrices, N.J. Higham, SIAM, 2008.
  /// \param A tensor
  /// \return \f$ \exp A \f$
  ///
  template<typename T>
  Tensor<T>
  exp_pade(Tensor<T> const & A);

  ///
  /// Logarithmic map by Taylor series, converges for \f$ |A-I| < 1 \f$
  /// \param A tensor
  /// \return \f$ \log A \f$
  ///
  template<typename T>
  Tensor<T>
  log(Tensor<T> const & A);

  ///
  /// Logarithmic map of a rotation
  /// \param R with \f$ R \in SO(3) \f$
  /// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
  ///
  template<typename T>
  Tensor<T>
  log_rotation(Tensor<T> const & R);

  ///
  /// Logarithmic map of a 180 degree rotation
  /// \param R with \f$ R \in SO(3) \f$
  /// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
  ///
  template<typename T>
  Tensor<T>
  log_rotation_pi(Tensor<T> const & R);

  /// Gaussian Elimination with partial pivot
  /// \param A
  /// \return \f$ xvec \f$
  ///
  template<typename T>
  Tensor<T>
  gaussian_elimination(Tensor<T> const & A);

  /// Apply Givens-Jacobi rotation on the left in place.
  /// \param c and s for a rotation G in form [c, s; -s, c]
  /// \param s
  /// \param i
  /// \param k
  /// \param A
  ///
  template<typename T>
  void
  givens_left(T const & c, T const & s, Index i, Index k, Tensor<T> & A);

  /// Apply Givens-Jacobi rotation on the right in place.
  /// \param A
  /// \param c and s for a rotation G in form [c, s; -s, c]
  /// \param s
  /// \param i
  /// \param k
  ///
  template<typename T>
  void
  givens_right(T const & c, T const & s, Index i, Index k, Tensor<T> & A);

  ///
  /// Exponential map of a skew-symmetric tensor
  /// \param r \f$ r \in so(3) \f$
  /// \return \f$ R = \exp R \f$ with \f$ R \in SO(3) \f$
  ///
  template<typename T>
  Tensor<T>
  exp_skew_symmetric(Tensor<T> const & r);

  ///
  /// Off-diagonal norm. Useful for SVD and other algorithms
  /// that rely on Jacobi-type procedures.
  /// \param A
  /// \return \f$ \sqrt(\sum_i \sum_{j, j\neq i} a_{ij}^2) \f$
  ///
  template<typename T>
  T
  norm_off_diagonal(Tensor<T> const & A);

  ///
  /// Arg max abs. Useful for inverse and other algorithms
  /// that rely on Jacobi-type procedures.
  /// \param A
  /// \return \f$ (p,q) = arg max_{i,j} |a_{ij}| \f$
  ///
  template<typename T>
  std::pair<Index, Index>
  arg_max_abs(Tensor<T> const & A);

  ///
  /// Arg max off-diagonal. Useful for SVD and other algorithms
  /// that rely on Jacobi-type procedures.
  /// \param A
  /// \return \f$ (p,q) = arg max_{i \neq j} |a_{ij}| \f$
  ///
  template<typename T>
  std::pair<Index, Index>
  arg_max_off_diagonal(Tensor<T> const & A);

  ///
  /// Sort and index. Useful for ordering singular values
  /// and eigenvalues and corresponding vectors in the
  /// respective decompositions.
  /// \param u vector to sort
  /// \return pair<v, P>
  /// \return v sorted vector
  /// \return P permutation matrix such that v = P^T u
  ///
  template<typename T>
  std::pair<Vector<T>, Tensor<T> >
  sort_permutation(Vector<T> const & u);

  ///
  /// Singular value decomposition (SVD)
  /// \param A tensor
  /// \return \f$ A = USV^T\f$
  ///
  template<typename T>
  boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
  svd(Tensor<T> const & A);

  ///
  /// Project to O(N) (Orthogonal Group) using a Newton-type algorithm.
  /// See Higham's Functions of Matrices p210 [2008]
  /// \param A tensor (often a deformation-gradient-like tensor)
  /// \return \f$ R = \arg min_Q \|A - Q\|\f$
  /// This algorithm projects a given tensor in GL(N) to O(N).
  /// The rotation/reflection obtained through this projection is
  /// the orthogonal component of the real polar decomposition
  ///
  template<typename T>
  Tensor<T>
  polar_rotation(Tensor<T> const & A);

  ///
  /// Left polar decomposition
  /// \param A tensor (often a deformation-gradient-like tensor)
  /// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
  ///
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_left(Tensor<T> const & A);

  ///
  /// Right polar decomposition
  /// \param A tensor (often a deformation-gradient-like tensor)
  /// \return \f$ RU = A \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
  ///
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_right(Tensor<T> const & A);

  ///
  /// Left polar decomposition computed with eigenvalue decomposition
  /// \param A tensor (often a deformation-gradient-like tensor)
  /// \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
  ///
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_left_eig(Tensor<T> const & A);

  ///
  /// R^3 right polar decomposition
  /// \param A tensor (often a deformation-gradient-like tensor)
  /// \return \f$ RU = F \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
  ///
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_right_eig(Tensor<T> const & A);

  ///
  /// Left polar decomposition with matrix logarithm for V
  /// \param F tensor (often a deformation-gradient-like tensor)
  /// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
  ///
  template<typename T>
  boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
  polar_left_logV(Tensor<T> const & F);

  ///
  /// Left polar decomposition with matrix logarithm for V using eig_spd_cos
  /// \param F tensor (often a deformation-gradient-like tensor)
  /// \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD, and log V
  ///
  template<typename T>
  boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
  polar_left_logV_lame(Tensor<T> const & F);

  ///
  /// Logarithmic map using BCH expansion (4 terms)
  /// \param v tensor
  /// \param r tensor
  /// \return Baker-Campbell-Hausdorff series up to 4 terms
  ///
  template<typename T>
  Tensor<T>
  bch(Tensor<T> const & v, Tensor<T> const & r);

  ///
  /// Symmetric Schur algorithm for R^2.
  /// \param \f$ A = [f, g; g, h] \in S(2) \f$
  /// \param f
  /// \param g
  /// \param h
  /// \return \f$ c, s \rightarrow [c, -s; s, c]\f$ diagonalizes A$
  ///
  template<typename T>
  std::pair<T, T>
  schur_sym(const T f, const T g, const T h);

  ///
  /// Givens rotation. [c, -s; s, c] [a; b] = [r; 0]
  /// \param a, b
  /// \return c, s
  ///
  template<typename T>
  std::pair<T, T>
  givens(T const & a, T const & b);

  ///
  /// Eigenvalue decomposition for symmetric 2nd-order tensor
  /// \param A tensor
  /// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
  ///
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  eig_sym(Tensor<T> const & A);

  ///
  /// Eigenvalue decomposition for SPD 2nd-order tensor
  /// \param A tensor
  /// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
  ///
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  eig_spd(Tensor<T> const & A);

  ///
  /// Eigenvalue decomposition for SPD 2nd-order tensor
  /// \param A tensor
  /// \return V eigenvectors, D eigenvalues in diagonal Matlab-style
  /// This algorithm comes from the journal article
  /// Scherzinger and Dohrmann, CMAME 197 (2008) 4007-4015
  ///
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  eig_spd_cos(Tensor<T> const & A);

  ///
  /// Cholesky decomposition, rank-1 update algorithm
  /// (Matrix Computations 3rd ed., Golub & Van Loan, p145)
  /// \param A assumed symetric tensor
  /// \return G Cholesky factor A = GG^T
  /// \return completed (bool) algorithm ran to completion
  ///
  template<typename T>
  std::pair<Tensor<T>, bool >
  cholesky(Tensor<T> const & A);

} // namespace Intrepid

#include "Intrepid_MiniTensor_LinearAlgebra.i.cc"
#include "Intrepid_MiniTensor_LinearAlgebra.t.cc"

#endif // Intrepid_MiniTensor_LinearAlgebra_h
