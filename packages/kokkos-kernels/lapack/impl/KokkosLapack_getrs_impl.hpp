// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_IMPL_GETRS_HPP_
#define KOKKOSLAPACK_IMPL_GETRS_HPP_

/// \file KokkosLapack_getrs_impl.hpp
/// \brief Implementation(s) of solve using LU factors.

#include <KokkosKernels_config.h>
#include <KokkosKernels_ArithTraits.hpp>
#include <KokkosBlas3_trsm.hpp>

namespace KokkosLapack {
namespace Impl {

// Tags for the laswp_functor forward/backward modes
struct ForwardMode {};
struct BackwardMode {};

// Functor parallelizing over columns:
//  - For square matrices this is probably fine
//  - For multi-vector it will slowdown with few vectors
template <class IpivView, class Matrix>
struct laswp_functor {
  Matrix m_A;
  IpivView m_Ipiv;

  laswp_functor(const Matrix& A, const IpivView& Ipiv) : m_A(A), m_Ipiv(Ipiv) {}

  void KOKKOS_FUNCTION operator()(ForwardMode, const int colIdx) const {
    typename Matrix::non_const_value_type tmp;
    for (int rowIdx = 0; rowIdx < m_A.extent_int(0); ++rowIdx) {
      const int piv       = m_Ipiv(rowIdx) - 1;
      tmp                 = m_A(rowIdx, colIdx);
      m_A(rowIdx, colIdx) = m_A(piv, colIdx);
      m_A(piv, colIdx)    = tmp;
    }
  }

  void KOKKOS_FUNCTION operator()(BackwardMode, const int colIdx) const {
    typename Matrix::non_const_value_type tmp;
    for (int rowIdx = m_A.extent_int(0) - 1; rowIdx >= 0; --rowIdx) {
      const int piv       = m_Ipiv(rowIdx) - 1;
      tmp                 = m_A(rowIdx, colIdx);
      m_A(rowIdx, colIdx) = m_A(piv, colIdx);
      m_A(piv, colIdx)    = tmp;
    }
  }
};

// Hinting at a potenatial future laswp function
// Right now not the actual focus of this PR
template <class ExecutionSpace, class IpivView, class Matrix>
void laswp_impl(const ExecutionSpace& space, const Matrix& A, const IpivView& Ipiv, const bool forward) {
  laswp_functor swapper(A, Ipiv);
  if (forward) {
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace, ForwardMode>(space, 0, A.extent(1)), swapper);
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace, BackwardMode>(space, 0, A.extent(1)), swapper);
  }
}

template <class ExecutionSpace, class AMatrix, class IpivView, class BMatrix, class InfoView>
void getrs_impl(const ExecutionSpace& space, const char trans[], const AMatrix& A, const IpivView& Ipiv,
                const BMatrix& B, const InfoView& /* Info */) {
  auto one = KokkosKernels::ArithTraits<typename AMatrix::non_const_value_type>::one();

  if (trans[0] == 'N' || trans[0] == 'n') {
    laswp_impl(space, B, Ipiv, true);
    KokkosBlas::trsm(space, "L", "L", "N", "U", one, A, B);
    KokkosBlas::trsm(space, "L", "U", "N", "N", one, A, B);
  } else {
    KokkosBlas::trsm(space, "L", "U", trans, "N", one, A, B);
    KokkosBlas::trsm(space, "L", "L", trans, "U", one, A, B);
    laswp_impl(space, B, Ipiv, false);
  }
}

}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_IMPL_GETRS_HPP_
