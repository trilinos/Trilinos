// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS2_SERIAL_GEMV_HPP_
#define KOKKOSBLAS2_SERIAL_GEMV_HPP_

#include "KokkosBlas2_serial_gemv_impl.hpp"
#include "KokkosBlas_util.hpp"

namespace KokkosBlas {
namespace Experimental {

template <class AlgoTag, class MatrixType, class XVector, class YVector, class ScalarType>
void KOKKOS_INLINE_FUNCTION serial_gemv(const char trans, const ScalarType& alpha, const MatrixType& A,
                                        const XVector& x, const ScalarType& beta, const YVector& y) {
  if (trans == 'N' || trans == 'n') {
    using mode = KokkosBlas::Trans::NoTranspose;
    KokkosBlas::SerialGemv<mode, AlgoTag>::invoke(alpha, A, x, beta, y);
  } else if (trans == 'T' || trans == 't') {
    using mode = KokkosBlas::Trans::Transpose;
    KokkosBlas::SerialGemv<mode, AlgoTag>::invoke(alpha, A, x, beta, y);
  } else if (trans == 'C' || trans == 'c') {
    using mode = KokkosBlas::Trans::ConjTranspose;
    KokkosBlas::SerialGemv<mode, AlgoTag>::invoke(alpha, A, x, beta, y);
  } else {
    Kokkos::abort("Matrix mode not supported");
  }
}

// default AlgoTag
template <class MatrixType, class XVector, class YVector, class ScalarType>
void KOKKOS_INLINE_FUNCTION serial_gemv(const char trans, const ScalarType& alpha, const MatrixType& A,
                                        const XVector& x, const ScalarType& beta, const YVector& y) {
  serial_gemv<KokkosBlas::Algo::Gemv::Default>(trans, alpha, A, x, beta, y);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
