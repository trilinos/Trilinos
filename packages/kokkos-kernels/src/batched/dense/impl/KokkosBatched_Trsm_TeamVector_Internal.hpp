#ifndef __KOKKOSBATCHED_TRSM_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSM_TEAMVECTOR_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

namespace KokkosBatched {

  ///
  /// Team Internal Impl
  /// ====================

  template<typename AlgoType>
  struct TeamVectorTrsmInternalLeftLower {
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const bool use_unit_diag,
           const int m, const int n, 
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ B, const int bs0, const int bs1);
  };

  template<>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  TeamVectorTrsmInternalLeftLower<Algo::Trsm::Unblocked>::
  invoke(const MemberType &member, 
         const bool use_unit_diag,
         const int m, const int n,
         const ScalarType alpha,
         const ValueType *__restrict__ A, const int as0, const int as1,
         /**/  ValueType *__restrict__ B, const int bs0, const int bs1) {

    const ScalarType one(1.0), zero(0.0);

    if (alpha == zero)   TeamVectorSetInternal  ::invoke(member, m, n, zero,  B, bs0, bs1);
    else {
      if (alpha != one)  TeamVectorScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
      if (m <= 0 || n <= 0) return 0;

      for (int p=0;p<m;++p) {
        // Made this non-const in order to WORKAROUND issue #349
        int iend = m-p-1;
        int jend = n;
          
        const ValueType
          *__restrict__ a21 = iend ? A+(p+1)*as0+p*as1 : NULL;
            
        ValueType
          *__restrict__ b1t =        B+p*bs0,
          *__restrict__ B2  = iend ? B+(p+1)*bs0 : NULL;

        member.team_barrier();
        if (!use_unit_diag) {
          const ValueType alpha11 = A[p*as0+p*as1];
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member,0,jend),[&](const int &j) {
              b1t[j*bs1] = b1t[j*bs1] / alpha11;
            });
          member.team_barrier();
        }
        Kokkos::parallel_for
	  (Kokkos::TeamThreadRange(member,iend),
	   [&](const int &i) {
	     Kokkos::parallel_for
	       (Kokkos::ThreadVectorRange(member,jend),
		[&](const int &j) {
		  // assume layout right for batched computation
		  B2[i*bs0+j*bs1] -= a21[i*as0] * b1t[j*bs1];
		});
	   });
      }
    }      
    return 0;
  }

  template<typename AlgoType>
  struct TeamVectorTrsmInternalLeftUpper {
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const bool use_unit_diag,
           const int m, const int n, 
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ B, const int bs0, const int bs1);
  };

  template<>
  template<typename MemberType,
           typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  TeamVectorTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::
  invoke(const MemberType &member, 
         const bool use_unit_diag,
         const int m, const int n,
         const ScalarType alpha,
         const ValueType *__restrict__ A, const int as0, const int as1,
         /**/  ValueType *__restrict__ B, const int bs0, const int bs1) {

    const ScalarType one(1.0), zero(0.0);

    // note that parallel range is different ( m*n vs m-1*n);        
    if (alpha == zero)  TeamVectorSetInternal  ::invoke(member, m, n, zero,  B, bs0, bs1);
    else {
      if (alpha != one) TeamVectorScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
      if (m <= 0 || n <= 0) return 0;
        
      ValueType *__restrict__ B0 = B;
      for (int p=(m-1);p>=0;--p) {
        // Made this non-const in order to WORKAROUND issue #349
        int iend = p;
        int jend = n;

        const ValueType *__restrict__ a01 = A+p*as1;
        /**/  ValueType *__restrict__ b1t = B+p*bs0;
            
        member.team_barrier();
        if (!use_unit_diag) {
          const ValueType alpha11 = A[p*as0+p*as1];
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member,0,jend),[&](const int &j) {
              b1t[j*bs1] = b1t[j*bs1] / alpha11;
            });
          member.team_barrier();
        }

        Kokkos::parallel_for
	  (Kokkos::TeamThreadRange(member,iend),
	   [&](const int &i) {
	     Kokkos::parallel_for
	       (Kokkos::ThreadVectorRange(member,jend),
		[&](const int &j) {
		  B0[i*bs0+j*bs1] -= a01[i*as0] * b1t[j*bs1];
		});
	   });
      }
    }
    return 0;
  }


}


#endif
