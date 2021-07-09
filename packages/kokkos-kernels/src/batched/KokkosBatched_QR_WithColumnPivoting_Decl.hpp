#ifndef __KOKKOSBATCHED_QR_WITH_COLUMNPIVOTING_DECL_HPP__
#define __KOKKOSBATCHED_QR_WITH_COLUMNPIVOTING_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// TeamVector QR
  ///

  template<typename MemberType,
           typename ArgAlgo>
  struct TeamVectorQR_WithColumnPivoting {
    template<typename AViewType,
             typename tViewType,
	     typename pViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const tViewType &t,
	   const pViewType &p,
           const wViewType &w,
	   /* */ int &matrix_rank);
  };

}

#include "KokkosBatched_QR_WithColumnPivoting_TeamVector_Impl.hpp"

#endif
