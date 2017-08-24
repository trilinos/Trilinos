#ifndef __KOKKOSBATCHED_TRSV_DECL_HPP__
#define __KOKKOSBATCHED_TRSV_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)


namespace KokkosBatched {
  namespace Experimental {  
    ///
    /// Serial Trsv
    ///

    template<typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct SerialTrsv {

      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b);
    };

    ///
    /// Team Trsv
    ///

    template<typename MemberType, 
             typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct TeamTrsv {

      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b);
    };

  }
}

#endif
