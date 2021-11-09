/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS3_TRMM_IMPL_HPP_
#define KOKKOSBLAS3_TRMM_IMPL_HPP_

/**
 * \file KokkosBlas3_trmm_impl.hpp
 * \brief Implementation of triangular matrix multiply
 */

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"
#include "KokkosBatched_Trmm_Decl.hpp"
#include "KokkosBatched_Trmm_Serial_Impl.hpp"

namespace KokkosBlas {
namespace Impl {

template <class AViewType, class BViewType>
void SerialTrmm_Invoke(const char side[], const char uplo[], const char trans[],
                       const char /*diag*/[],
                       typename BViewType::const_value_type& alpha,
                       const AViewType& A, const BViewType& B) {
  using KokkosBatched::Algo;
  using KokkosBatched::Diag;
  using KokkosBatched::SerialTrmmInternalLeftLower;
  using KokkosBatched::SerialTrmmInternalLeftUpper;
  using KokkosBatched::SerialTrmmInternalRightLower;
  using KokkosBatched::SerialTrmmInternalRightUpper;

  char __side = tolower(side[0]), __uplo = tolower(uplo[0]),
       __trans = tolower(trans[0]);
  //__diag = tolower(diag[0]);
  bool do_conj = true;

  // Ignoring diag, see "ech-note" in KokkosBatched_Trmm_Serial_Internal.hpp

  //// Lower non-transpose ////
  if (__side == 'l' && __uplo == 'l' && __trans == 'n')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1),
        B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'l' && __trans == 'n')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1),
        B.data(), B.stride(0), B.stride(1));
  //// Lower transpose /////
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 't')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'l' && __trans == 't')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));

  //// Lower conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 'c')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'l' && __trans == 'c')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));
  //// Upper non-transpose ////
  if (__side == 'l' && __uplo == 'u' && __trans == 'n')
    SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1),
        B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'u' && __trans == 'n')
    SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(0), A.extent(1),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1),
        B.data(), B.stride(0), B.stride(1));
  //// Upper transpose
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 't')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'u' && __trans == 't')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, !do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));

  //// Upper conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 'c')
    SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));
  if (__side == 'r' && __uplo == 'u' && __trans == 'c')
    SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, do_conj, A.extent(1), A.extent(0),
        B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0),
        B.data(), B.stride(0), B.stride(1));
}
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSBLAS3_TRMM_IMPL_HPP_
