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

#ifndef KOKKOSCOMPAT_VIEW_DEF_HPP
#define KOKKOSCOMPAT_VIEW_DEF_HPP

#include "KokkosCompat_View.hpp"

namespace Kokkos {
  namespace Compat {

    // template <typename D, typename T>
    // Kokkos::View<T*,D>
    // getKokkosViewDeepCopy(const Teuchos::ArrayView<T>& a) {
    //   typedef typename std::conditional<
    //         Impl::VerifyExecutionCanAccessMemorySpace< D, Kokkos::HostSpace>::value,
    //         typename D::execution_space, Kokkos::HostSpace>::type
    //       HostDevice;
    //   typedef Kokkos::View<T*,D>  view_type;
    //   typedef Kokkos::View<T*,typename view_type::array_layout,HostDevice,Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
    //   if (a.size() == 0)
    //     return view_type();
    //   view_type v("", a.size());
    //   unmanaged_host_view_type hv(a.getRawPtr(), a.size());
    //   Kokkos::deep_copy(v,hv);
    //   return v;
    // }

    // template <typename D, typename T>
    // Kokkos::View<const T*,D>
    // getKokkosViewDeepCopy(const Teuchos::ArrayView<const T>& a) {
    //   typedef typename std::conditional<
    //         Impl::VerifyExecutionCanAccessMemorySpace< D, Kokkos::HostSpace>::value,
    //         typename D::execution_space, Kokkos::HostSpace>::type
    //       HostDevice;
    //   typedef Kokkos::View<T*,D>  view_type;
    //   typedef Kokkos::View<const T*,typename view_type::array_layout,HostDevice,Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
    //   if (a.size() == 0)
    //     return view_type();
    //   view_type v("", a.size());
    //   unmanaged_host_view_type hv(a.getRawPtr(), a.size());
    //   Kokkos::deep_copy(v,hv);
    //   return v;
    // }

    // template <typename T, typename D>
    // void
    // realloc(Kokkos::View<T*,D>& v, const typename D::size_type size) {
    //   Kokkos::realloc(v,size);
    // }

    // template <typename T, typename L, typename D, typename M>
    // Kokkos::View<T*,L,D,M>
    // create_view(const std::string& label, size_t size) {
    //   return Kokkos::View<T*,L,D,M>(label, size);
    // }

    // #define COMPAT_INSTANT(T,D)

    // template Kokkos::View<T*,D> getKokkosViewDeepCopy<D,T>(const Teuchos::ArrayView<T>& a);
    //template Kokkos::View<const T*,D> getKokkosViewDeepCopy<D,T>(const Teuchos::ArrayView<const T>& a);
    // template void realloc<T,D>(Kokkos::View<T*,D>& v, const D::size_type size);
    // template Kokkos::View<T*,Kokkos::LayoutLeft,D,void> create_view<T,Kokkos::LayoutLeft,D,void>(const std::string& label, size_t size);
    // template Kokkos::View<T*,Kokkos::LayoutRight,D,void> create_view<T,Kokkos::LayoutRight,D,void>(const std::string& label, size_t size);
    // template Kokkos::View<T*,D,void,void> create_view<T,D,void,void>(const std::string& label, size_t size);

  } // namespace Compat
} // namespace Kokkos

#endif // KOKKOSCOMPAT_VIEW_DEF_HPP
