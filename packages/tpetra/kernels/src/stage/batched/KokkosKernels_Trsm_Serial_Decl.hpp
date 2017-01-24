#ifndef __KOKKOSKERNELS_TRSM_SERIAL_DECL_HPP__
#define __KOKKOSKERNELS_TRSM_SERIAL_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <immintrin.h>

namespace KokkosKernels {

  namespace Serial {

    // specialized for different m and n
    // Solve L(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTriSolveLeftLowerUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveLeftLowerUnitDiag(const int as0, const int as1,
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
    struct InnerTriSolveLeftLowerNonUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveLeftLowerNonUnitDiag(const int as0, const int as1,
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
    struct InnerTriSolveRightUpperUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveRightUpperUnitDiag(const int as0, const int as1,
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
    struct InnerTriSolveRightUpperNonUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveRightUpperNonUnitDiag(const int as0, const int as1,
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
             const AViewType A,
             /**/  BViewType B);
      
    };

    template<typename ArgSide,
             typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag>
    struct Trsm<ArgSide,ArgUplo,ArgTrans,ArgDiag,void> {

      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType A,
             /**/  BViewType B) {
        return 0;
      }
      
    };

  }

}

#endif
