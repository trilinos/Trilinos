// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_TRTRI_IMPL_HPP_
#define KOKKOSLAPACK_TRTRI_IMPL_HPP_

/**
 * \file KokkosLapack_trtri_impl.hpp
 * \brief Implementation of triangular matrix inverse
 */

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosBatched_Trtri_Decl.hpp"
#include "KokkosBatched_Trtri_Serial_Impl.hpp"

namespace KokkosLapack {
namespace Impl {

template <class RViewType, class AViewType>
void SerialTrtri_Invoke(const RViewType &R, const char uplo[], const char diag[], const AViewType &A) {
  using KokkosBatched::Algo;
  using KokkosBatched::Diag;
  using KokkosBatched::SerialTrtriInternalLower;
  using KokkosBatched::SerialTrtriInternalUpper;

  char _uplo = tolower(uplo[0]), _diag = tolower(diag[0]);

  //// Lower ////
  if (_uplo == 'l') {
    if (_diag == 'u') {
      R() = SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(Diag::Unit::use_unit_diag, A.extent(0),
                                                                     A.extent(1), A.data(), A.stride(0), A.stride(1));
    } else {
      R() = SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(Diag::NonUnit::use_unit_diag, A.extent(0),
                                                                     A.extent(1), A.data(), A.stride(0), A.stride(1));
    }
  } else {
    //// Upper ////
    if (_diag == 'u') {
      R() = SerialTrtriInternalUpper<Algo::Trtri::Unblocked>::invoke(Diag::Unit::use_unit_diag, A.extent(0),
                                                                     A.extent(1), A.data(), A.stride(0), A.stride(1));
    } else {
      R() = SerialTrtriInternalUpper<Algo::Trtri::Unblocked>::invoke(Diag::NonUnit::use_unit_diag, A.extent(0),
                                                                     A.extent(1), A.data(), A.stride(0), A.stride(1));
    }
  }
}
}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSLAPACK_TRTRI_IMPL_HPP_
