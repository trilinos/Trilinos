#ifndef __KOKKOSBATCHED_COPY_DECL_HPP__
#define __KOKKOSBATCHED_COPY_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)


namespace KokkosBatched {
  namespace Experimental {
    ///
    /// Serial Copy
    ///

    template<typename ArgTrans>
    struct SerialCopy {
      template<typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const AViewType &A,
             /* */ BViewType &B);
    };

    ///
    /// Team Copy
    ///

    template<typename MemberType, typename ArgTrans>
    struct TeamCopy {
      template<typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const AViewType &A,
             const BViewType &B);
    };

  }
}


#endif
