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
#ifndef KOKKOSKERNELS_HELPERS_HPP_
#define KOKKOSKERNELS_HELPERS_HPP_

#include "KokkosKernels_config.h"  // KOKKOSKERNELS_INST_LAYOUTLEFT, KOKKOSKERNELS_INST_LAYOUTRIGHT
#include "KokkosKernels_default_types.hpp"  // default_layout

namespace KokkosKernels {
namespace Impl {

// Unify Layout of a View to PreferredLayoutType if possible
// (either matches already, or is rank-0/rank-1 and contiguous)
// Used to reduce number of code instantiations.
template <class ViewType, class PreferredLayoutType>
struct GetUnifiedLayoutPreferring {
  typedef typename std::conditional<
      ((ViewType::rank == 1) && (!std::is_same<typename ViewType::array_layout,
                                               Kokkos::LayoutStride>::value)) ||
          ((ViewType::rank == 0)),
      PreferredLayoutType, typename ViewType::array_layout>::type array_layout;
};

template <class ViewType>
struct GetUnifiedLayout {
  using array_layout =
      typename GetUnifiedLayoutPreferring<ViewType,
                                          default_layout>::array_layout;
};

template <class T, class TX, bool do_const,
          bool isView = Kokkos::is_view<T>::value>
struct GetUnifiedScalarViewType {
  typedef typename TX::non_const_value_type type;
};

template <class T, class TX>
struct GetUnifiedScalarViewType<T, TX, false, true> {
  typedef Kokkos::View<typename T::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<
                           T, typename TX::array_layout>::array_layout,
                       typename T::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      type;
};

template <class T, class TX>
struct GetUnifiedScalarViewType<T, TX, true, true> {
  typedef Kokkos::View<typename T::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<
                           T, typename TX::array_layout>::array_layout,
                       typename T::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      type;
};

}  // namespace Impl
}  // namespace KokkosKernels
#endif
