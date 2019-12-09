// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_LOCAL_DEEP_COPY_HPP
#define TPETRA_DETAILS_LOCAL_DEEP_COPY_HPP

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include <type_traits>
#include <cstring> // memcpy

namespace Tpetra {
namespace Details {

  template<class DT, class ... DP , class ST, class ... SP>
  KOKKOS_INLINE_FUNCTION void
  local_deep_copy_contiguous (const Kokkos::View<DT, DP...>& dst,
                              const Kokkos::View<ST, SP...>& src)
  {
    auto dst_raw = dst.data ();
    auto src_raw = src.data ();

    using dst_type = Kokkos::View<DT, DP...>;
    using src_type = Kokkos::View<ST, SP...>;
    using dst_value_type = typename dst_type::non_const_value_type;
    using src_value_type = typename src_type::non_const_value_type;
    constexpr bool same_value =
      std::is_same<dst_value_type, src_value_type>::value;
    constexpr bool trivially_copyable =
      std::is_trivially_copyable<dst_value_type>::value;
    constexpr bool both_rank_one =
      static_cast<int> (dst_type::Rank) == 1 &&
      static_cast<int> (src_type::Rank) == 1;
    // All contiguous rank-1 Views have the "same layout."
    constexpr bool same_layout = both_rank_one ||
      std::is_same<typename dst_type::array_layout,
                   typename dst_type::array_layout>::value;
    if (same_value && trivially_copyable && same_layout) {
      // See Trilinos #6393; the casts avoid compiler warnings for the
      // "else" case.  Once we have "if constexpr" we can get rid of
      // the casts.
      std::memcpy (reinterpret_cast<char*> (dst_raw),
                   reinterpret_cast<const char*> (src_raw),
                   sizeof (dst_type));
    }
    else {
      for (size_t k = 0; k < dst.span (); ++k) {
        dst_raw[k] = src_raw[k];
      }
    }
  }

  template<class DT, class ... DP , class ST, class ... SP>
  KOKKOS_INLINE_FUNCTION void
  local_deep_copy (const Kokkos::View<DT, DP...>& dst,
                   const Kokkos::View<ST, SP...>& src,
                   typename std::enable_if<
                     static_cast<int> (Kokkos::View<DT, DP...>::Rank) == 1 &&
                     static_cast<int> (Kokkos::View<ST, SP...>::Rank) == 1
                   >::type* = 0)
  {
    // FIXME (mfh 08 Dec 2019) Use Kokkos::local_deep_copy, once
    // that function leaves the Kokkos::Experimental namespace.
    // using Kokkos::Experimental::local_deep_copy;
    // local_deep_copy (dst, src);

    if (dst.span_is_contiguous () && src.span_is_contiguous ()) {
      local_deep_copy_contiguous (dst, src);
    }
    else {
      using size_type = decltype (dst.extent (0));
      for (size_type k = 0; k < dst.extent (0); ++k) {
        dst(k) = src(k);
      }
    }
  }

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LOCAL_DEEP_COPY_HPP
