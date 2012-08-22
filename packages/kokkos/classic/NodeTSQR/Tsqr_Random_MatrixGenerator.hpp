//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

#ifndef __TSQR_Random_MatrixGenerator_hpp
#define __TSQR_Random_MatrixGenerator_hpp

#include <Tsqr_Blas.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_ScalarTraits.hpp>

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Random {

    template< class Ordinal, class Scalar, class Generator >
    class MatrixGenerator {
    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      typedef Generator generator_type;

      MatrixGenerator (Generator& generator) : gen_ (generator) {}

      void
      fill_random (const Ordinal nrows, 
		   const Ordinal ncols, 
		   Scalar A[], 
		   const Ordinal lda)
      {
	for (Ordinal j = 0; j < ncols; ++j)
	  {
	    Scalar* const A_j = &A[j*lda];
	    for (Ordinal i = 0; i < nrows; ++i)
	      A_j[i] = gen_();
	  }
      }

      /// Fill the nrows by ncols matrix Q (in column-major order, with
      /// leading dimension ldq >= nrows) with a random orthogonal
      /// matrix.
      ///
      /// \note If you want the resulting Q factor to be from a Haar
      /// distribution (the usual one for random orthogonal matrices),
      /// Generator must have a normal distribution.
      void
      explicit_Q (const Ordinal nrows, 
		  const Ordinal ncols, 
		  Scalar Q[], 
		  const Ordinal ldq)
      {
	// Fill Q with random numbers
	this->fill_random (nrows, ncols, Q, ldq);

	// Get ready for QR factorization
	LAPACK< Ordinal, Scalar > lapack;
	std::vector< Scalar > tau (std::min(nrows, ncols));

	// Workspace query
	Scalar _lwork1, _lwork2;
	int info = 0;
	lapack.GEQRF (nrows, ncols, Q, ldq, &tau[0], &_lwork1, -1, &info);
	if (info != 0)
	  throw std::logic_error("LAPACK GEQRF LWORK query failed");

	lapack.ORGQR (nrows, ncols, ncols, Q, ldq, &tau[0], &_lwork2, -1, &info);
	if (info != 0)
	  throw std::logic_error("LAPACK ORGQR LWORK query failed");

	// Allocate workspace.  abs() returns a magnitude_type, and we
	// can compare those using std::max.  If Scalar is complex,
	// you can't compare it using max.
	const Ordinal lwork = checkedCast (std::max (ScalarTraits< Scalar >::abs (_lwork1), 
						     ScalarTraits< Scalar >::abs (_lwork2)));
	std::vector< Scalar > work (lwork);

	// Factor the input matrix
	lapack.GEQRF (nrows, ncols, Q, ldq, &tau[0], &work[0], lwork, &info);
	if (info != 0)
	  throw std::runtime_error("LAPACK GEQRF failed");

	// Compute explicit Q factor in place
	lapack.ORGQR (nrows, ncols, ncols, Q, ldq, &tau[0], &work[0], lwork, &info);
	if (info != 0)
	  throw std::runtime_error("LAPACK ORGQR failed");
      }


      /// Fill the nrows by ncols matrix Q (in column-major order, with
      /// leading dimension ldq >= nrows) with a random orthogonal
      /// matrix, stored implicitly.  tau (of length min(nrows,ncols))
      /// is part of this storage.
      ///
      /// \note If you want the resulting Q factor to be from a Haar
      /// distribution (the usual one for random orthogonal matrices),
      /// Generator must have a normal distribution.
      void
      implicit_Q (const Ordinal nrows, 
		  const Ordinal ncols, 
		  Scalar Q[], 
		  const Ordinal ldq,
		  Scalar tau[])
      {
	// Fill Q with random numbers
	this->fill_random (nrows, ncols, Q, ldq);

	// Get ready for QR factorization
	LAPACK< Ordinal, Scalar > lapack;

	// Workspace query
	Scalar _lwork1;
	int info = 0;
	lapack.GEQRF (nrows, ncols, Q, ldq, tau, &_lwork1, -1, &info);
	if (info != 0)
	  throw std::logic_error("LAPACK GEQRF LWORK query failed");

	// Allocate workspace.
	const Ordinal lwork = checkedCast (ScalarTraits< Scalar >::abs (_lwork1));
	std::vector< Scalar > work (lwork);

	// Factor the input matrix
	lapack.GEQRF (nrows, ncols, Q, ldq, tau, &work[0], lwork, &info);
	if (info != 0)
	  throw std::runtime_error("LAPACK GEQRF failed");
      }

      template< class MatrixViewType >
      void
      implicit_Q (MatrixViewType& Q, 
		  typename MatrixViewType::scalar_type tau[])
      {
	implicit_Q (Q.nrows(), Q.ncols(), Q.get(), Q.lda(), tau);
      }

      void
      fill_random_svd (const Ordinal nrows, 
		       const Ordinal ncols, 
		       Scalar A[], 
		       const Ordinal lda,
		       const magnitude_type singular_values[])
      {
	typedef Matrix< Ordinal, Scalar > matrix_type;
	typedef MatView< Ordinal, Scalar > matrix_view_type;

	matrix_type U (nrows, ncols, Scalar(0));
	matrix_type V (ncols, ncols, Scalar(0));
	std::vector<Scalar> tau_U (std::min (nrows, ncols));
	std::vector<Scalar> tau_V (ncols);

	// Fill A with zeros, and then make its diagonal the given set
	// of singular values.
	matrix_view_type A_view (nrows, ncols, A, lda);
	A_view.fill (Scalar (0));
	for (Ordinal j = 0; j < ncols; ++j)
	  // Promote magnitude_type to Scalar here.
	  A_view(j,j) = Scalar (singular_values[j]);

	// Generate random orthogonal U (nrows by ncols) and V (ncols by
	// ncols).  Keep them stored implicitly.
	implicit_Q (U, &tau_U[0]);
	implicit_Q (V, &tau_V[0]);

	// Workspace query for ORMQR.
	Scalar _lwork1, _lwork2;
	int info = 0;
	LAPACK< Ordinal, Scalar > lapack;
	lapack.ORMQR ("L", "N", nrows, ncols, ncols, U.get(), U.lda(), &tau_U[0], 
		      A, lda, &_lwork1, -1, &info);
	if (info != 0)
	  {
	    std::ostringstream os;
	    os << "LAPACK ORMQR LWORK query failed with INFO = " << info 
	       << ": called ORMQR(\"L\", \"N\", " << nrows << ", " << ncols 
	       << ", " << ncols << ", NULL, " << U.lda() << ", NULL, NULL, " 
	       << lda << ", WORK, -1, &INFO)";
	    throw std::logic_error(os.str());
	  }
	if (ScalarTraits< Scalar >::is_complex)
	  lapack.ORMQR ("R", "C", nrows, ncols, ncols, V.get(), V.lda(), &tau_V[0], 
			A, lda, &_lwork2, -1, &info);
	else
	  lapack.ORMQR ("R", "T", nrows, ncols, ncols, V.get(), V.lda(), &tau_V[0], 
			A, lda, &_lwork2, -1, &info);
	if (info != 0)
	  throw std::logic_error("LAPACK ORMQR LWORK query failed");

	// Allocate workspace.
	Ordinal lwork = checkedCast (std::max (ScalarTraits< Scalar >::abs (_lwork1), 
					       ScalarTraits< Scalar >::abs (_lwork2)));
	std::vector< Scalar > work (lwork);

	// Apply U to the left side of A, and V^H to the right side of A.
	lapack.ORMQR ("L", "N", nrows, ncols, ncols, U.get(), U.lda(), &tau_U[0], 
		      A, lda, &work[0], lwork, &info);
	if (info != 0)
	  throw std::runtime_error("LAPACK ORMQR failed (first time)");
	if (ScalarTraits< Scalar >::is_complex)
	  lapack.ORMQR ("R", "C", nrows, ncols, ncols, V.get(), V.lda(), &tau_V[0], 
			A, lda, &work[0], lwork, &info);
	else
	  lapack.ORMQR ("R", "T", nrows, ncols, ncols, V.get(), V.lda(), &tau_V[0], 
			A, lda, &work[0], lwork, &info);
	if (info != 0)
	  throw std::runtime_error("LAPACK ORMQR failed (second time)");
      }


      /// \brief Fill in a random upper triangular matrix
      ///
      /// Fill R with a random n by n upper triangular matrix, with
      /// the specified singular values.
      ///
      /// \param n [in] Dimension of R (n by n)
      ///
      /// \param R [out] On output: n by n upper triangular matrix
      ///   with stride ldr, stored in non-packed, column-major
      ///   format.
      ///
      /// \param ldr [in] Stride (a.k.a. leading dimension) of R
      ///
      /// \param singular_values [in] n nonnegative real values, to be
      ///   used as the singular values of the output R.
      void
      fill_random_R (const Ordinal n,
		     Scalar R[],
		     const Ordinal ldr,
		     const magnitude_type singular_values[])
      {
	// Fill R with an n x n (not upper triangular) random matrix
	// having the given singular values.
	fill_random_svd (n, n, R, ldr, singular_values);

	// Compute the QR factorization in place of R (which isn't upper triangular yet).
	std::vector< Scalar > tau (n);

	// Workspace size query for QR factorization.
	Scalar _lwork1;
	int info = 0;
	LAPACK< Ordinal, Scalar > lapack;
	lapack.GEQRF (n, n, R, ldr, &tau[0], &_lwork1, -1, &info);
	if (info != 0)
	  throw std::logic_error("LAPACK GEQRF LWORK query failed");

	// Allocate workspace
	Ordinal lwork = checkedCast (ScalarTraits< Scalar >::abs (_lwork1));
	std::vector< Scalar > work (lwork);

	// Compute QR factorization (implicit representation in place).
	lapack.GEQRF (n, n, R, ldr, &tau[0], &work[0], lwork, &info);
	if (info != 0)
	  throw std::runtime_error("LAPACK GEQRF failed");

	// Zero out the stuff below the diagonal of R, leaving just the R factor.
	for (Ordinal j = 0; j < n; ++j)
	  for (Ordinal i = j+1; i < n; ++i)
	    R[i + j*ldr] = Scalar(0);
      }

    private:
      static Ordinal 
      checkedCast (const magnitude_type& x)
      {
	if (x < std::numeric_limits< Ordinal >::min() || x > std::numeric_limits< Ordinal >::max())
	  throw std::range_error("Scalar input cannot be safely cast to an Ordinal");
	else if (std::numeric_limits< magnitude_type >::is_signed && 
		 x < magnitude_type(0) &&
		 ! std::numeric_limits< Ordinal >::is_signed)
	  throw std::range_error("Scalar input is negative, but Ordinal is unsigned");
	else
	  return static_cast< Ordinal > (x);
      }

      Generator& gen_;
    };
  } // namespace Random
} // namespace TSQR

#endif // __TSQR_Random_MatrixGenerator_hpp
