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

#ifndef KOKKOSBLAS_TRTRI_IMPL_HPP_
#define KOKKOSBLAS_TRTRI_IMPL_HPP_

/**
 * \file KokkosBlas_trtri_impl.hpp
 * \brief Implementation of triangular matrix inverse
 */


#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosBatched_Trtri_Decl.hpp"
#include "KokkosBatched_Trtri_Serial_Impl.hpp"

using namespace KokkosBatched;

namespace KokkosBlas {
  namespace Impl {

    template<class RViewType,
             class AViewType>
    void SerialTrtri_Invoke (const RViewType &R,
                            const char uplo[],
                            const char diag[],
                            const AViewType &A)
    { 
      char __uplo = tolower(uplo[0]),
           __diag = tolower(diag[0]);

      //// Lower ////
      if (__uplo == 'l') {
        if (__diag == 'u') {
          R() = SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                                  A.extent(0), A.extent(1),
                                                                  A.data(), A.stride(0), A.stride(1));
        } else {
          R() = SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                                  A.extent(0), A.extent(1),
                                                                  A.data(), A.stride(0), A.stride(1));
        }
      } else {
      //// Upper ////
        if (__diag == 'u') {
          R() = SerialTrtriInternalUpper<Algo::Trtri::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                                  A.extent(0), A.extent(1),
                                                                  A.data(), A.stride(0), A.stride(1));
        } else {
          R() = SerialTrtriInternalUpper<Algo::Trtri::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                                  A.extent(0), A.extent(1),
                                                                  A.data(), A.stride(0), A.stride(1));
        }
      }
    }
  } // namespace Impl
} // namespace KokkosBlas
#endif // KOKKOSBLAS_TRTRI_IMPL_HPP_
