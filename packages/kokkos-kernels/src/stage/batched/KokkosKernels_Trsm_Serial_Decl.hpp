#ifndef __KOKKOSKERNELS_TRSM_SERIAL_DECL_HPP__
#define __KOKKOSKERNELS_TRSM_SERIAL_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <immintrin.h>

namespace KokkosKernels {

  namespace Serial {

    // specialized for different m and n
    // Solve L(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsmLeftLowerUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTrsmLeftLowerUnitDiag(const int as0, const int as1,
                                 const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int n,
                 /**/  ValueType *__restrict__ B);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };

    // specialized for different m and n
    // Solve L(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsmLeftLowerNonUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTrsmLeftLowerNonUnitDiag(const int as0, const int as1,
                                    const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int n,
                 /**/  ValueType *__restrict__ B);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };

    // specialized for different m and n
    // Solve  X(m x n) U(n x n) = B(m x n)
    template<int bmn>
    struct InnerTrsmRightUpperUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTrsmRightUpperUnitDiag(const int as0, const int as1,
                                  const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m,
                 /**/  ValueType *__restrict__ B);
      
      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };    

    // specialized for different m and n
    // Solve  X(m x n) U(n x n) = B(m x n)
    template<int bmn>
    struct InnerTrsmRightUpperNonUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTrsmRightUpperNonUnitDiag(const int as0, const int as1,
                                     const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m,
                 /**/  ValueType *__restrict__ B);
      
      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };    

    // specialized for different m and n
    // Solve U(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsmLeftUpperUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTrsmLeftUpperUnitDiag(const int as0, const int as1,
                                 const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int n,
                 /**/  ValueType *__restrict__ B);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };

    // specialized for different m and n
    // Solve U(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTrsmLeftUpperNonUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTrsmLeftUpperNonUnitDiag(const int as0, const int as1,
                                    const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int n,
                 /**/  ValueType *__restrict__ B);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };

    template<typename ArgSide,
             typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct Trsm {

      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        static_assert(false, "KokkosKernels::Trsm::invoke:: Not yet implemented");
        return 0;
      }
      
    };

  }

}

#endif
