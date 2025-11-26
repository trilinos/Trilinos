// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS3_TRMM_IMPL_HPP_
#define KOKKOSBLAS3_TRMM_IMPL_HPP_

/**
 * \file KokkosBlas3_trmm_impl.hpp
 * \brief Implementation of triangular matrix multiply
 */

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosKernels_ArithTraits.hpp"
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

  char _side = tolower(side[0]), _uplo = tolower(uplo[0]), _trans = tolower(trans[0]);
  //__diag = tolower(diag[0]);
  bool do_conj = true;

  // Ignoring diag, see "ech-note" in KokkosBatched_Trmm_Serial_Internal.hpp

  //// Lower non-transpose ////
  if (_side == 'l' && _uplo == 'l' && _trans == 'n')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  if (_side == 'r' && _uplo == 'l' && _trans == 'n')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  //// Lower transpose /////
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (_side == 'l' && _uplo == 'l' && _trans == 't')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (_side == 'r' && _uplo == 'l' && _trans == 't')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));

  //// Lower conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (_side == 'l' && _uplo == 'l' && _trans == 'c')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (_side == 'r' && _uplo == 'l' && _trans == 'c')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  //// Upper non-transpose ////
  if (_side == 'l' && _uplo == 'u' && _trans == 'n')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  if (_side == 'r' && _uplo == 'u' && _trans == 'n')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
  //// Upper transpose
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (_side == 'l' && _uplo == 'u' && _trans == 't')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (_side == 'r' && _uplo == 'u' && _trans == 't')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));

  //// Upper conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (_side == 'l' && _uplo == 'u' && _trans == 'c')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (_side == 'r' && _uplo == 'u' && _trans == 'c')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
}
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSBLAS3_TRMM_IMPL_HPP_
