#ifndef __KOKKOSKERNELS_LU_SERIAL_DECL_HPP__
#define __KOKKOSKERNELS_LU_SERIAL_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosKernels {

  namespace Serial {

    template<int bmn>
    struct InnerLU {
      const int _as0, _as1;

      KOKKOS_INLINE_FUNCTION
      InnerLU(const int as0, const int as1) 
        : _as0(as0), _as1(as1) {}

      // lu
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(ValueType *__restrict__ A);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const int m, const int n,
                 ValueType *__restrict__ A);
    };


    template<typename ArgAlgo>
    struct LU {

      // no piv version
      template<typename AViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const AViewType &A) {
        static_assert(false, "KokkosKernels::LU:: Not yet implemented");
        return 0;
      }

    };

  }

}

#endif
