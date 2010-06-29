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

  template< class LocalOrdinal, class Scalar >
  Scalar
  global_frobenius_norm (const LocalOrdinal nrows_local, 
			 const LocalOrdinal ncols,
			 const Scalar A_local[],
			 const LocalOrdinal lda_local,
			 MessengerBase< Scalar >* messenger)
  {
    // FIXME (mfh 20 Apr 2010) This is currently implemented using an
    // all-reduction.  This may result in different processors getting
    // slightly different answers, due to floating-point arithmetic
    // roundoff.  We might not want this if we are using this function
    // to test a routine.

    Scalar local_result (0);
    for (LocalOrdinal j = 0; j < ncols; j++)
      {
	const Scalar* const cur_col = &A_local[j*lda_local];
	for (LocalOrdinal i = 0; i < nrows_local; ++i)
	  {
	    const Scalar xi = cur_col[i];
	    local_result = local_result + xi * xi;
	  }
      }
    Scalar global_result = messenger->globalSum (local_result);
    return sqrt (global_result);
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
		 MessengerBase< Scalar >* messenger)
  {
    using std::make_pair;
    using std::pair;
    using std::vector;
  
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

    const int my_rank = messenger->rank();
    BLAS< LocalOrdinal, Scalar > blas;

    //
    // Compute $\| I - Q^T * Q \|_F$
    //

    // Compute Q_local^T * Q_local (this node's component of Q^T*Q)
    vector< Scalar > Temp (ncols*ncols, std::numeric_limits< Scalar >::quiet_NaN());
    const LocalOrdinal ld_temp = ncols;
    blas.GEMM ("T", "N", ncols, ncols, nrows_local, 
	       1.0, Q_local, ldq_local, Q_local, ldq_local, 
	       0.0, &Temp[0], ld_temp);
  
    // Reduce over all the processors to get the global Q^T*Q in Temp2.
    vector< Scalar > Temp2 (ncols*ncols, std::numeric_limits< Scalar >::quiet_NaN());
    messenger->globalVectorSum (&Temp[0], &Temp2[0], ncols*ncols);

    // Compute I-(Q^T*Q) redundantly on all processors
    for (LocalOrdinal j = 0; j < ncols; j++)
      Temp2[j + j*ld_temp] = 1.0 - Temp2[j + j*ld_temp];
  
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
	       -1.0, Q_local, ldq_local, R, ldr, 
	       1.0, &Resid[0], ld_resid);

    const magnitude_type Resid_F = 
      global_frobenius_norm (nrows_local, ncols, &Resid[0], ld_resid, messenger);

    return make_pair (Resid_F / A_F, Orthog_F / A_F);
  }

} // namespace TSQR

#endif // __TSQR_Tsqr_GlobalVerify_hpp

