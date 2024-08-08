//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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

  char __uplo = tolower(uplo[0]), __diag = tolower(diag[0]);

  //// Lower ////
  if (__uplo == 'l') {
    if (__diag == 'u') {
      R() = SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(Diag::Unit::use_unit_diag, A.extent(0),
                                                                     A.extent(1), A.data(), A.stride(0), A.stride(1));
    } else {
      R() = SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(Diag::NonUnit::use_unit_diag, A.extent(0),
                                                                     A.extent(1), A.data(), A.stride(0), A.stride(1));
    }
  } else {
    //// Upper ////
    if (__diag == 'u') {
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
