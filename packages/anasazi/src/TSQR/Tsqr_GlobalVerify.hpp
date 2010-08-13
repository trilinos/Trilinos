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

#ifndef __TSQR_Tsqr_GlobalVerify_hpp
#define __TSQR_Tsqr_GlobalVerify_hpp

#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Blas.hpp>
#include <Tsqr_Util.hpp>

#include <utility> // std::pair
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // Unfortunately, you need c++0x support to have default template
  // arguments of template functions.  Otherwise we would make this a
  // template function and set the default value of isComplex to
  // ScalarTraits< Scalar >::is_complex.  Also, C++ doesn't like
  // partial specialization of template functions, for no good reason.
  // So we have to make this a class.
  template< class Scalar, bool isComplex = ScalarTraits< Scalar >::is_complex >
  class GlobalSummer {
  public:
    typedef Scalar scalar_type;
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

    static magnitude_type
    sum (const Scalar& localSum,
	 MessengerBase< Scalar >* const messenger);
  };

  // Complex-arithmetic "forward declaration"
  template< class Scalar >
  class GlobalSummer< Scalar, true > {
  public:
    typedef Scalar scalar_type;
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

    static magnitude_type
    sum (const Scalar& localSum,
	 MessengerBase< Scalar >* const messenger);
  };

  // Real-arithmetic "forward declaration"
  template< class Scalar >
  class GlobalSummer< Scalar, false > {
  public:
    typedef Scalar scalar_type;
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

    static magnitude_type
    sum (const Scalar& localSum,
	 MessengerBase< Scalar >* const messenger);
  };

  // Complex-arithmetic case
  template< class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  GlobalSummer< Scalar, true >::sum (const Scalar& localSum,
				     MessengerBase< Scalar >* const messenger)
  {
    // In order to use a MessengerBase<Scalar> on magnitude_type
    // values, we have to convert local_result to a Scalar, and then
    // convert back the result.  We convert by setting the real
    // component of the Scalar to the magnitude_type.  This isn't
    // guaranteed to work if magnitude_type has a greater dynamic
    // range than Scalar.  That's possible, but that's not how we do
    // things with ScalarTraits< std::complex< T > >, and that's not
    // how LAPACK does it either, so it's fair to assume that
    // magnitude_type and the individual components of Scalar have the
    // same dynamic range.
    const magnitude_type localSumAbs = ScalarTraits< Scalar >::abs (localSum);
    const Scalar localSumAsScalar (localSumAbs, magnitude_type(0));
    const Scalar globalSumAsScalar = messenger->globalSum (localSumAsScalar);
    const magnitude_type globalSum = ScalarTraits< Scalar >::abs (globalSumAsScalar);
    return globalSum;
  }

  // Real-arithmetic case
  template< class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  GlobalSummer< Scalar, false >::sum (const Scalar& localSum,
				      MessengerBase< Scalar >* const messenger)
  {
    const Scalar localSumAsScalar (localSum);
    const Scalar globalSumAsScalar = messenger->globalSum (localSumAsScalar);
    const magnitude_type globalSum = ScalarTraits< Scalar >::abs (globalSumAsScalar);
    return globalSum;
  }

  template< class LocalOrdinal, class Scalar >
  typename ScalarTraits< Scalar >::magnitude_type
  global_frobenius_norm (const LocalOrdinal nrows_local, 
			 const LocalOrdinal ncols,
			 const Scalar A_local[],
			 const LocalOrdinal lda_local,
			 MessengerBase< Scalar >* const messenger)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

    // FIXME (mfh 20 Apr 2010) This is currently implemented using an
    // all-reduction.  This may result in different processors getting
    // slightly different answers, due to floating-point arithmetic
    // roundoff.  We might not want this if we are using this function
    // to test a routine.

    magnitude_type local_result (0);
    for (LocalOrdinal j = 0; j < ncols; j++)
      {
	const Scalar* const cur_col = &A_local[j*lda_local];
	for (LocalOrdinal i = 0; i < nrows_local; ++i)
	  {
	    const magnitude_type abs_xi = ScalarTraits< Scalar >::abs (cur_col[i]);
	    local_result = local_result + abs_xi * abs_xi;
	  }
      }
    // GlobalSummmer() is a hack to let us use a Scalar - type
    // MessengerBase with magnitude_type inputs and outputs.
    // Otherwise we would need to carry around a MessengerBase<
    // magnitude_type > object as well.
    const magnitude_type globalResult = 
      GlobalSummer< Scalar, ScalarTraits< Scalar >::is_complex >::sum (local_result, messenger);
    return sqrt (globalResult);
  }

  template< class LocalOrdinal, class Scalar >
  std::pair< typename ScalarTraits< Scalar >::magnitude_type, typename ScalarTraits< Scalar >::magnitude_type >
  global_verify (const LocalOrdinal nrows_local, 
		 const LocalOrdinal ncols, 
		 const Scalar A_local[],
		 const LocalOrdinal lda_local,
		 const Scalar Q_local[],
		 const LocalOrdinal ldq_local,
		 const Scalar R[],
		 const LocalOrdinal ldr,
		 MessengerBase< Scalar >* const messenger)
  {
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    using std::make_pair;
    using std::pair;
    using std::vector;

    const magnitude_type ZERO (0);
    const magnitude_type ONE (1);
    BLAS< LocalOrdinal, Scalar > blas;

    //
    // Compute $\| I - Q^T * Q \|_F$
    //

    // Compute Q_local^T * Q_local (this node's component of Q^T*Q)
    vector< Scalar > Temp (ncols*ncols, std::numeric_limits< Scalar >::quiet_NaN());
    const LocalOrdinal ld_temp = ncols;

    if (ScalarTraits< Scalar >::is_complex)
      blas.GEMM ("C", "N", ncols, ncols, nrows_local, 
		 ONE, Q_local, ldq_local, Q_local, ldq_local, 
		 ZERO, &Temp[0], ld_temp);
    else
      blas.GEMM ("T", "N", ncols, ncols, nrows_local, 
		 ONE, Q_local, ldq_local, Q_local, ldq_local, 
		 ZERO, &Temp[0], ld_temp);
  
    // Reduce over all the processors to get the global Q^T*Q in Temp2.
    vector< Scalar > Temp2 (ncols*ncols, std::numeric_limits< Scalar >::quiet_NaN());
    messenger->globalVectorSum (&Temp[0], &Temp2[0], ncols*ncols);

    // Compute I-(Q^T*Q) redundantly on all processors
    for (LocalOrdinal j = 0; j < ncols; j++)
      Temp2[j + j*ld_temp] = ONE - Temp2[j + j*ld_temp];
  
    // Compute the Frobenius norm of I - Q^T*Q, redundantly on all processors.
    const magnitude_type Orthog_F = 
      local_frobenius_norm (ncols, ncols, &Temp2[0], ld_temp);

    // Compute the Frobenius norm of A.  
    const magnitude_type A_F = 
      global_frobenius_norm (nrows_local, ncols, &A_local[0], lda_local, messenger);

    //
    // Compute $\| A - Q*R \|_F$
    //

    vector< Scalar > Resid (nrows_local * ncols, std::numeric_limits< Scalar >::quiet_NaN());
    const LocalOrdinal ld_resid = nrows_local;

    // Resid := A (deep copy)
    copy_matrix (nrows_local, ncols, &Resid[0], ld_resid, A_local, lda_local);

    // Resid := Resid - Q*R
    blas.GEMM ("N", "N", nrows_local, ncols, ncols, 
	       -ONE, Q_local, ldq_local, R, ldr, 
	       ONE, &Resid[0], ld_resid);

    const magnitude_type Resid_F = 
      global_frobenius_norm (nrows_local, ncols, &Resid[0], ld_resid, messenger);

    return make_pair (Resid_F / A_F, Orthog_F / A_F);
  }

} // namespace TSQR

#endif // __TSQR_Tsqr_GlobalVerify_hpp

