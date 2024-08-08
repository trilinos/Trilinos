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

#ifndef KOKKOSBLAS3_TRMM_IMPL_HPP_
#define KOKKOSBLAS3_TRMM_IMPL_HPP_

/**
 * \file KokkosBlas3_trmm_impl.hpp
 * \brief Implementation of triangular matrix multiply
 */

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosBatched_Trmm_Decl.hpp"
#include "KokkosBatched_Trmm_Serial_Impl.hpp"

namespace KokkosBlas {
namespace Impl {

template <class AViewType, class BViewType>
void SerialTrmm_Invoke(const char side[], const char uplo[], const char trans[], const char /*diag*/[],
                       typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B) {
  using KokkosBatched::Algo;
  using KokkosBatched::Diag;
  using KokkosBatched::SerialTrmmInternalLeftLower;
  using KokkosBatched::SerialTrmmInternalLeftUpper;
  using KokkosBatched::SerialTrmmInternalRightLower;
  using KokkosBatched::SerialTrmmInternalRightUpper;

  char __side = tolower(side[0]), __uplo = tolower(uplo[0]), __trans = tolower(trans[0]);
  //__diag = tolower(diag[0]);
  bool do_conj = true;

  // Ignoring diag, see "ech-note" in KokkosBatched_Trmm_Serial_Internal.hpp

  //// Lower non-transpose ////
  if (__side == 'l' && __uplo == 'l' && __trans == 'n')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'l' && __trans == 'n')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  //// Lower transpose /////
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 't')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'l' && __trans == 't')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));

  //// Lower conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 'c')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'l' && __trans == 'c')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  //// Upper non-transpose ////
  if (__side == 'l' && __uplo == 'u' && __trans == 'n')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'u' && __trans == 'n')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  //// Upper transpose
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 't')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'u' && __trans == 't')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));

  //// Upper conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 'c')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'u' && __trans == 'c')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
}
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSBLAS3_TRMM_IMPL_HPP_
