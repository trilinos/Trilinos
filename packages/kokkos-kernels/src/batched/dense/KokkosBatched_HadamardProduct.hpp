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
#ifndef __KOKKOSBATCHED_HADAMARDPRODUCT_HPP__
#define __KOKKOSBATCHED_HADAMARDPRODUCT_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// \brief Serial Batched Hadamard Product:
///   v_ij <- x_ij * y_ij for all i = 1, ..., n and j = 1, ..., N
/// where:
///   * n is the number of rows,
///   * N is the number of vectors.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam VViewType: Input type for V, needs to be a 2D view
///
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in]: Input vector Y, a rank 2 view
/// \param V [out]: Output vector V, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///

struct SerialHadamardProduct {
  template <typename XViewType, typename YViewType, typename VViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &X,
                                           const YViewType &Y,
                                           const VViewType &V);
};

/// \brief Team Batched Hadamard Product:
///   v_ij <- x_ij * y_ij for all i = 1, ..., n and j = 1, ..., N
/// where:
///   * n is the number of rows,
///   * N is the number of vectors.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam VViewType: Input type for V, needs to be a 2D view
///
/// \param member [in]: TeamPolicy member
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in]: Input vector Y, a rank 2 view
/// \param V [out]: Output vector V, a rank 2 view
///
/// A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType>
struct TeamHadamardProduct {
  template <typename XViewType, typename YViewType, typename VViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const XViewType &X,
                                           const YViewType &Y,
                                           const VViewType &V);
};

/// \brief TeamVector Batched Hadamard Product:
///   v_ij <- x_ij * y_ij for all i = 1, ..., n and j = 1, ..., N
/// where:
///   * n is the number of rows,
///   * N is the number of vectors.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam VViewType: Input type for V, needs to be a 2D view
///
/// \param member [in]: TeamPolicy member
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in]: Input vector Y, a rank 2 view
/// \param V [out]: Output vector V, a rank 2 view
///
/// Two nested parallel_for with both TeamThreadRange and ThreadVectorRange
/// (or one with TeamVectorRange) are used inside.
///

template <typename MemberType>
struct TeamVectorHadamardProduct {
  template <typename XViewType, typename YViewType, typename VViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const XViewType &X,
                                           const YViewType &Y,
                                           const VViewType &V);
};

template <typename MemberType, typename ArgMode>
struct HadamardProduct {
  template <typename XViewType, typename YViewType, typename VViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const XViewType &X,
                                           const YViewType &Y,
                                           const VViewType &V) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialHadamardProduct::template invoke<XViewType, YViewType,
                                                     VViewType>(X, Y, V);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val =
          TeamHadamardProduct<MemberType>::template invoke<XViewType, YViewType,
                                                           VViewType>(member, X,
                                                                      Y, V);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      r_val = TeamVectorHadamardProduct<MemberType>::template invoke<
          XViewType, YViewType, VViewType>(member, X, Y, V);
    }
    return r_val;
  }
};
}  // namespace KokkosBatched

#include "KokkosBatched_HadamardProduct_Impl.hpp"

#endif
