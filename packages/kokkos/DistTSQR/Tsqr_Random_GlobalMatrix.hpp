// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef __Tsqr_Random_GlobalMatrix_hpp
#define __Tsqr_Random_GlobalMatrix_hpp

#include "Tsqr_Blas.hpp"
#include "Tsqr_Matrix.hpp"
#include "Tsqr_Random_MatrixGenerator.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tsqr_RMessenger.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Random {

    template<class MatrixViewType>
    static void
    scaleMatrix (MatrixViewType& A, 
		 const typename MatrixViewType::scalar_type& denom)
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      typedef typename MatrixViewType::scalar_type scalar_type;

      const ordinal_type nrows = A.nrows();
      const ordinal_type ncols = A.ncols();
      const ordinal_type lda = A.lda();

      if (nrows == lda)
	{ // The whole A matrix is stored contiguously.
	  const ordinal_type nelts = nrows * ncols;
	  scalar_type* const A_ptr = A.get();
	  std::transform (A_ptr, A_ptr + nelts, A_ptr, std::bind2nd(std::divides<scalar_type>(), denom));
	}
      else
	{ // Each column of A is stored contiguously.
	  for (ordinal_type j = 0; j < ncols; ++j)
	    {
	      scalar_type* const A_j = &A(0,j);
	      std::transform (A_j, A_j + nrows, A_j, std::bind2nd(std::divides<scalar_type>(), denom));
	    }
	}
    }

    template< class MatrixViewType, class Generator >
    void
    randomGlobalMatrix (Generator* const pGenerator, 
			MatrixViewType& A_local,
			const typename Teuchos::ScalarTraits< typename MatrixViewType::scalar_type >::magnitudeType singular_values[],
			MessengerBase< typename MatrixViewType::ordinal_type >* const ordinalMessenger,
			MessengerBase< typename MatrixViewType::scalar_type >* const scalarMessenger)
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      typedef typename MatrixViewType::scalar_type scalar_type;
      typedef typename Teuchos::ScalarTraits< scalar_type >::magnitudeType magnitude_type;
      using std::vector;

      const bool b_local_debug = false;

      const int rootProc = 0;
      const int nprocs = ordinalMessenger->size();
      const int myRank = ordinalMessenger->rank();
      BLAS< ordinal_type, scalar_type > blas;

      const ordinal_type nrowsLocal = A_local.nrows();
      const ordinal_type ncols = A_local.ncols();

      // Theory: Suppose there are P processors.  Proc q wants an m_q by n
      // component of the matrix A, which we write as A_q.  On Proc 0, we
      // generate random m_q by n orthogonal matrices Q_q (in explicit
      // form), and send Q_q to Proc q.  The m by n matrix [Q_0; Q_1; ...;
      // Q_{P-1}] is not itself orthogonal.  However, the m by n matrix
      // Q = [Q_0 / P; Q_1 / P; ...; Q_{P-1} / P] is orthogonal:
      // 
      // \sum_{q = 0}^{P-1} (Q_q^T * Q_q) / P = I.

      if (myRank == rootProc)
	{
	  typedef Random::MatrixGenerator< ordinal_type, scalar_type, Generator > matgen_type;
	  matgen_type matGen (*pGenerator);

	  // Generate a random ncols by ncols upper triangular matrix
	  // R with the given singular values.
	  Matrix< ordinal_type, scalar_type > R (ncols, ncols, scalar_type(0));
	  matGen.fill_random_R (ncols, R.get(), R.lda(), singular_values);

	  // Broadcast R to all the processors.
	  scalarMessenger->broadcast (R.get(), ncols*ncols, rootProc);

	  // Generate (for myself) a random nrowsLocal x ncols
	  // orthogonal matrix, stored in explicit form.
	  Matrix< ordinal_type, scalar_type > Q_local (nrowsLocal, ncols);
	  matGen.explicit_Q (nrowsLocal, ncols, Q_local.get(), Q_local.lda());

	  // Scale the (local) orthogonal matrix by the number of
	  // processors P, to make the columns of the global matrix Q
	  // orthogonal.  (Otherwise the norm of each column will be P
	  // instead of 1.)
	  const scalar_type P = static_cast< scalar_type > (nprocs);
	  // Do overflow check.  If casting P back to scalar_type
	  // doesn't produce the same value as nprocs, the cast
	  // overflowed.  We take the real part, because scalar_type
	  // might be complex.
	  if (nprocs != static_cast<int> (Teuchos::ScalarTraits<scalar_type>::real (P)))
	    throw std::runtime_error ("Casting nprocs to Scalar failed");

	  scaleMatrix (Q_local, P);

	  // A_local := Q_local * R
	  blas.GEMM ("N", "N", nrowsLocal, ncols, ncols,
		     scalar_type(1), Q_local.get(), Q_local.lda(), 
		     R.get(), R.lda(), 
		     scalar_type(0), A_local.get(), A_local.lda());

	  for (int recvProc = 1; recvProc < nprocs; ++recvProc)
	    {
	      // Ask the receiving processor how big (i.e., how many rows)
	      // its local component of the matrix is.
	      ordinal_type nrowsRemote = 0;
	      ordinalMessenger->recv (&nrowsRemote, 1, recvProc, 0);

	      if (b_local_debug)
		{
		  std::ostringstream os;
		  os << "For Proc " << recvProc << ": local block is " 
		     << nrowsRemote << " by " << ncols << std::endl;
		  std::cerr << os.str();
		}

	      // Make sure Q_local is big enough to hold the data for
	      // the current receiver proc.
	      Q_local.reshape (nrowsRemote, ncols);

	      // Compute a random nrowsRemote * ncols orthogonal
	      // matrix Q_local, for the current receiving processor.
	      matGen.explicit_Q (nrowsRemote, ncols, Q_local.get(), Q_local.lda());

	      // Send Q_local to the current receiving processor.
	      scalarMessenger->send (Q_local.get(), nrowsRemote*ncols, recvProc, 0);
	    }
	}
      else
	{
	  // Receive the R factor from Proc 0.  There's only 1 R
	  // factor for all the processes.
	  Matrix< ordinal_type, scalar_type > R (ncols, ncols, scalar_type (0));
	  scalarMessenger->broadcast (R.get(), ncols*ncols, rootProc);

	  // Q_local (nrows_local by ncols, random orthogonal matrix)
	  // will be received from Proc 0, where it was generated.
	  const ordinal_type recvSize = nrowsLocal * ncols;
	  Matrix< ordinal_type, scalar_type > Q_local (nrowsLocal, ncols);

	  // Tell Proc 0 how many rows there are in the random orthogonal
	  // matrix I want to receive from Proc 0.
	  ordinalMessenger->send (&nrowsLocal, 1, rootProc, 0);

	  // Receive the orthogonal matrix from Proc 0.
	  scalarMessenger->recv (Q_local.get(), recvSize, rootProc, 0);

	  // Scale the (local) orthogonal matrix by the number of
	  // processors, to make the global matrix Q orthogonal.
	  const scalar_type P = static_cast< scalar_type > (nprocs);
	  // Do overflow check.  If casting P back to scalar_type
	  // doesn't produce the same value as nprocs, the cast
	  // overflowed.  We take the real part, because scalar_type
	  // might be complex.
	  if (nprocs != static_cast<int> (Teuchos::ScalarTraits<scalar_type>::real (P)))
	    throw std::runtime_error ("Casting nprocs to Scalar failed");
	  scaleMatrix (Q_local, P);

	  // A_local := Q_local * R
	  blas.GEMM ("N", "N", nrowsLocal, ncols, ncols, 
		     scalar_type(1), Q_local.get(), Q_local.lda(),
		     R.get(), R.lda(),
		     scalar_type(0), A_local.get(), A_local.lda());
	}
    }
  } // namespace Random
} // namespace TSQR

#endif // __Tsqr_Random_GlobalMatrix_hpp
