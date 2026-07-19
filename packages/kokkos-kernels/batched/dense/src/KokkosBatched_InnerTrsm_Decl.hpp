// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_INNER_TRSM_DECL_HPP
#define KOKKOSBATCHED_INNER_TRSM_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

// specialized for different m and n
// Solve L(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftLowerUnitDiag {
  const int m_as0, m_as1, m_bs0, m_bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftLowerUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : m_as0(as0), m_as1(as1), m_bs0(bs0), m_bs1(bs1) {}

  // trisolve
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

// specialized for different m and n
// Solve L(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftLowerNonUnitDiag {
  const int m_as0, m_as1, m_bs0, m_bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftLowerNonUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : m_as0(as0), m_as1(as1), m_bs0(bs0), m_bs1(bs1) {}

  // trisolve
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

// specialized for different m and n
// Solve U(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftUpperUnitDiag {
  const int m_as0, m_as1, m_bs0, m_bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftUpperUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : m_as0(as0), m_as1(as1), m_bs0(bs0), m_bs1(bs1) {}

  // trisolve
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

// specialized for different m and n
// Solve U(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftUpperNonUnitDiag {
  const int m_as0, m_as1, m_bs0, m_bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftUpperNonUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : m_as0(as0), m_as1(as1), m_bs0(bs0), m_bs1(bs1) {}

  // trisolve
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename Op, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(Op op, const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

}  // namespace KokkosBatched

#endif
