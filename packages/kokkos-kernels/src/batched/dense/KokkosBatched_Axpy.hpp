//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.4
//       Copyright (2021) National Technology & Engineering
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
#ifndef __KOKKOSBATCHED_AXPY_HPP__
#define __KOKKOSBATCHED_AXPY_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// \brief Serial Batched AXPY:
///   y_l <- alpha_l * x_l + y_l for all l = 1, ..., N
/// where:
///   * N is the number of vectors,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in/out]: Output vector Y, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///

struct SerialAxpy {
  template <typename XViewType, typename YViewType, typename alphaViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const alphaViewType &alpha,
                                           const XViewType &X,
                                           const YViewType &Y);
};

/// \brief Team Batched AXPY:
///   y_l <- alpha_l * x_l + y_l for all l = 1, ..., N
/// where:
///   * N is the number of vectors,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param member [in]: TeamPolicy member
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in/out]: Output vector Y, a rank 2 view
///
/// A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType>
struct TeamAxpy {
  template <typename XViewType, typename YViewType, typename alphaViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const alphaViewType &alpha,
                                           const XViewType &X,
                                           const YViewType &Y);
};

/// \brief TeamVector Batched AXPY:
///   y_l <- alpha_l * x_l + y_l for all l = 1, ..., N
/// where:
///   * N is the number of vectors,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param member [in]: TeamPolicy member
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in/out]: Output vector Y, a rank 2 view
///
/// Two nested parallel_for with both TeamThreadRange and ThreadVectorRange
/// (or one with TeamVectorRange) are used inside.
///

template <typename MemberType>
struct TeamVectorAxpy {
  template <typename XViewType, typename YViewType, typename alphaViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const alphaViewType &alpha,
                                           const XViewType &X,
                                           const YViewType &Y);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Axpy_Impl.hpp"

#endif
