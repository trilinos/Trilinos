#ifndef __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Gemv_TeamVector_Internal.hpp"
#include "KokkosBatched_Trsv_TeamVector_Internal.hpp"

#include "KokkosBatched_Gemm_TeamVector_Internal.hpp"
#include "KokkosBatched_Trsm_TeamVector_Internal.hpp"

namespace KokkosBatched {

  ///
  /// TeamVector Internal
  /// =================== 
  struct TeamVectorSolveUTV_Internal {
    template<typename MemberType,
             typename ValueType,
	     typename IntType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
	   const int matrix_rank,
           const int m,
           const ValueType * U, const int us0, const int us1,
	   const ValueType * T, const int ts0, const int ts1,
	   const ValueType * V, const int vs0, const int vs1,
	   const IntType   * p, const int ps0,
	   /* */ ValueType * x, const int xs0,
	   /* */ ValueType * b, const int bs0,
           /* */ ValueType * w) {
      typedef ValueType value_type;
      //typedef IntType int_type;

      const value_type one(1), zero(0);
      const int ws0 = 1;

      if (matrix_rank < m) {
	/// w = U^T b
	TeamVectorGemvInternal<Algo::Gemv::Unblocked>
	  ::invoke(member,
		   matrix_rank, m,
		   one,
		   U, us1, us0,
		   b, bs0,
		   zero,
		   w, ws0);
	
	/// w = T^{-1} w
	TeamVectorTrsvInternalLower<Algo::Trsv::Unblocked>
	  ::invoke(member,
		   false,
		   matrix_rank,
		   one,
		   T, ts0, ts1,
		   w, ws0);
	
	/// x = V^T w
	TeamVectorGemvInternal<Algo::Gemv::Unblocked>
	  ::invoke(member,
		   m, matrix_rank, 
		   one,
		   V, vs1, vs0,
		   w, ws0,
		   zero,
		   x, xs0);
      } else {
	TeamVectorGemvInternal<Algo::Gemv::Unblocked>
	  ::invoke(member,
		   matrix_rank, m,
		   one,
		   U, us1, us0,
		   b, bs0,
		   zero,
		   x, xs0);

	TeamVectorTrsvInternalUpper<Algo::Trsv::Unblocked>
	  ::invoke(member,
		   false,
		   matrix_rank,
		   one,
		   T, ts0, ts1,
		   x, xs0);	
      }

      /// x = P^T x
      TeamVectorApplyPivotVectorBackwardInternal
	::invoke(member,
		 m,
		 p, ps0,
		 x, xs0);

      return 0;
    }

    template<typename MemberType,
             typename ValueType,
	     typename IntType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
	   const int matrix_rank,
           const int m, const int nrhs,
           const ValueType * U, const int us0, const int us1,
	   const ValueType * T, const int ts0, const int ts1,
	   const ValueType * V, const int vs0, const int vs1,
	   const IntType   * p, const int ps0,
	   /* */ ValueType * X, const int xs0, const int xs1,
	   /* */ ValueType * B, const int bs0, const int bs1,
           /* */ ValueType * w) {
      typedef ValueType value_type;
      //typedef IntType int_type;

      const value_type one(1), zero(0);

      value_type * W = w; /// m x nrhs
      const int ws0 = xs0 < xs1 ? 1 : nrhs, ws1 = xs0 < xs1 ? m : 1;

      if (matrix_rank < m) {
	/// U is m x matrix_rank
	/// T is matrix_rank x matrix_rank
	/// V is matrix_rank m
	/// W = U^T B
	TeamVectorGemmInternal<Algo::Gemm::Unblocked>
	  ::invoke(member,
		   matrix_rank, nrhs, m,
		   one,
		   U, us1, us0,
		   B, bs0, bs1,
		   zero,
		   W, ws0, ws1);

	/// W = T^{-1} W
	TeamVectorTrsmInternalLeftLower<Algo::Trsm::Unblocked>
	  ::invoke(member,
		   false,
		   matrix_rank, nrhs,
		   one,
		   T, ts0, ts1,
		   W, ws0, ws1);
	
	/// X = V^T W
	TeamVectorGemmInternal<Algo::Gemm::Unblocked>
	  ::invoke(member,
		   m, nrhs, matrix_rank, 
		   one,
		   V, vs1, vs0,
		   W, ws0, ws1,
		   zero,
		   X, xs0, xs1);
      } else {
	TeamVectorGemmInternal<Algo::Gemm::Unblocked>
	  ::invoke(member,
		   m, nrhs, matrix_rank, 
		   one,
		   U, us1, us0,
		   B, bs0, bs1,
		   zero,
		   X, xs0, xs1);

	TeamVectorTrsmInternalLeftUpper<Algo::Trsm::Unblocked>
	  ::invoke(member,
		   false,
		   matrix_rank, nrhs,
		   one,
		   T, ts0, ts1,
		   X, xs0, xs1);		
      }
      
      /// X = P^T X
      TeamVectorApplyPivotMatrixBackwardInternal
      	::invoke(member,
      		 nrhs, m,
      		 p, ps0,
      		 X, xs0, xs1);

      return 0;
    }

  };

} // end namespace KokkosBatched


#endif
