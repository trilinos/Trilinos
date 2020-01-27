#ifndef __KOKKOSBATCHED_INNER_MULTIPLE_DOT_PRODUCT_DECL_HPP__
#define __KOKKOSBATCHED_INNER_MULTIPLE_DOT_PRODUCT_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {

  template<int mb>
  struct InnerMultipleDotProduct {
    const int _as0, _as1, _xs0, _ys0;

    KOKKOS_INLINE_FUNCTION
    InnerMultipleDotProduct(const int as0, const int as1,
                            const int xs0,
                            const int ys0) 
      : _as0(as0), _as1(as1),
        _xs0(xs0),
        _ys0(ys0) {}

    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int 
    serial_invoke(const ScalarType alpha,
                  const ValueType *__restrict__ A,
                  const ValueType *__restrict__ x,
                  const int n, 
                  /**/  ValueType *__restrict__ y);

    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    serial_invoke(const ScalarType alpha,
                  const ValueType *__restrict__ A,
                  const ValueType *__restrict__ x,
                  const int m, const int n, 
                  /**/  ValueType *__restrict__ y);
  };
}



#endif
