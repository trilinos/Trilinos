#ifndef __KOKKOSKERNELS_TRSV_DECL_HPP__
#define __KOKKOSKERNELS_TRSV_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <immintrin.h>

namespace KokkosKernels {

  namespace Serial {

    // specialized for different m and n
    // Solve L(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsvLowerUnitDiag {
      const int _as0, _as1, _bs0;

      KOKKOS_INLINE_FUNCTION
      InnerTrsvLowerUnitDiag(const int as0, const int as1,
                             const int bs0)
        : _as0(as0), _as1(as1),
          _bs0(bs0) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 /**/  ValueType *__restrict__ b);
      
      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, 
                 /**/  ValueType *__restrict__ b);
    };

    // specialized for different m and n
    // Solve L(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsvLowerNonUnitDiag {
      const int _as0, _as1, _bs0;
      
      KOKKOS_INLINE_FUNCTION
      InnerTrsvLowerNonUnitDiag(const int as0, const int as1,
                                const int bs0)
        : _as0(as0), _as1(as1),
          _bs0(bs0) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 /**/  ValueType *__restrict__ b);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, 
                 /**/  ValueType *__restrict__ b);
    };

    // specialized for different m and n
    // Solve U(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsvUpperUnitDiag {
      const int _as0, _as1, _bs0;

      KOKKOS_INLINE_FUNCTION
      InnerTrsvUpperUnitDiag(const int as0, const int as1,
                             const int bs0)
        : _as0(as0), _as1(as1),
          _bs0(bs0) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 /**/  ValueType *__restrict__ b);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, 
                 /**/  ValueType *__restrict__ b);
    };

    // specialized for different m and n
    // Solve U(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsvUpperNonUnitDiag {
      const int _as0, _as1, _bs0;
      
      KOKKOS_INLINE_FUNCTION
      InnerTrsvUpperNonUnitDiag(const int as0, const int as1,
                                const int bs0)
        : _as0(as0), _as1(as1),
          _bs0(bs0) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 /**/  ValueType *__restrict__ b);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, 
                 /**/  ValueType *__restrict__ b);
    };

    template<typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct Trsv {

      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        static_assert(false, "KokkosKernels::Trsv:: Not yet implemented");
        return 0;
      }
      
    };


  }

}

#endif
