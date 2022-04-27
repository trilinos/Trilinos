#ifndef __KOKKOSBATCHED_INNER_LU_DECL_HPP__
#define __KOKKOSBATCHED_INNER_LU_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {

template <int bmn>
struct InnerLU {
  const int _as0, _as1;

  KOKKOS_INLINE_FUNCTION
  InnerLU(const int as0, const int as1) : _as0(as0), _as1(as1) {}

  // lu
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(ValueType *KOKKOS_RESTRICT A);

  // for remainder square
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const int m,
                                           ValueType *KOKKOS_RESTRICT A);

  // for remainder
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const int m, const int n,
                                           ValueType *KOKKOS_RESTRICT A);
};
}  // namespace KokkosBatched

#endif
