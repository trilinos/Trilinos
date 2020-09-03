#ifndef __KOKKOSBATCHED_UTV_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_UTV_TEAMVECTOR_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_SetTriangular_Internal.hpp"
#include "KokkosBatched_QR_TeamVector_Internal.hpp"
#include "KokkosBatched_QR_WithColumnPivoting_TeamVector_Internal.hpp"
#include "KokkosBatched_QR_FormQ_TeamVector_Internal.hpp"

namespace KokkosBatched {

  ///
  /// TeamVector Internal
  /// =================== 
  struct TeamVectorUTV_Internal {
    template<typename MemberType,
             typename ValueType,
	     typename IntType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const int m, // m = NumRows(A)
           /* */ ValueType * A, const int as0, const int as1,
	   /* */ IntType   * p, const int ps0,
	   /* */ ValueType * U, const int us0, const int us1,
	   /* */ ValueType * V, const int vs0, const int vs1,
           /* */ ValueType * w, // 3*m, tau, norm, householder workspace
	   /* */ int &matrix_rank) {
      typedef ValueType value_type;
      //typedef IntType int_type;

      value_type *t = w; w+= m;
      const int ts0(1);

      value_type *work = w;

      matrix_rank = -1;
      TeamVectorQR_WithColumnPivotingInternal
      	::invoke(member,
      		 m, m,
      		 A, as0, as1,
      		 t, ts0,
      		 p, ps0,
      		 work,
      		 matrix_rank);
      
      TeamVectorQR_FormQ_Internal
      	::invoke(member,
      		 m, matrix_rank,
      		 A, as0, as1,
      		 t, ts0,
      		 U, us0, us1,
      		 work);

      /// for rank deficient matrix
      if (matrix_rank < m) {
	const value_type zero(0);
	TeamVectorSetLowerTriangularInternal
	  ::invoke(member,
		   matrix_rank, matrix_rank,
		   1, zero,
		   A, as0, as1);
	
	TeamVectorQR_Internal
	  ::invoke(member,
		   m, matrix_rank,
		   A, as1, as0,
		   t, ts0,
		   work);
	
	TeamVectorQR_FormQ_Internal
	  ::invoke(member,
		   m, matrix_rank,
		   A, as1, as0,
		   t, ts0,
		   V, vs1, vs0,
		   work);
      }

      return 0;
    }
  };

} // end namespace KokkosBatched


#endif
