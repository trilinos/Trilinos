/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSCOMPAT_VIEW_DEF_HPP
#define KOKKOSCOMPAT_VIEW_DEF_HPP

#include "KokkosCompat_View.hpp"

namespace Kokkos {
  namespace Compat {

    template <typename D, typename T>
    Kokkos::View<T*,D>
    getKokkosViewDeepCopy(const Teuchos::ArrayView<T>& a) {
      typedef typename Kokkos::Impl::if_c<
            VerifyExecutionSpaceCanAccessDataSpace< D, Kokkos::HostSpace>::value,
            typename D::execution_space, Kokkos::HostSpace>::type
          HostDevice;
      typedef Kokkos::View<T*,D>  view_type;
      typedef Kokkos::View<T*,typename view_type::array_layout,HostDevice,Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
      if (a.size() == 0)
        return view_type();
      view_type v("", a.size());
      unmanaged_host_view_type hv(a.getRawPtr(), a.size());
      Kokkos::deep_copy(v,hv);
      return v;
    }

    template <typename D, typename T>
    Kokkos::View<const T*,D>
    getKokkosViewDeepCopy(const Teuchos::ArrayView<const T>& a) {
      typedef typename Kokkos::Impl::if_c<
            VerifyExecutionSpaceCanAccessDataSpace< D, Kokkos::HostSpace>::value,
            typename D::execution_space, Kokkos::HostSpace>::type
          HostDevice;
      typedef Kokkos::View<T*,D>  view_type;
      typedef Kokkos::View<const T*,typename view_type::array_layout,HostDevice,Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
      if (a.size() == 0)
        return view_type();
      view_type v("", a.size());
      unmanaged_host_view_type hv(a.getRawPtr(), a.size());
      Kokkos::deep_copy(v,hv);
      return v;
    }

    template <typename T, typename D>
    void
    realloc(Kokkos::View<T*,D>& v, const typename D::size_type size) {
      Kokkos::realloc(v,size);
    }

    template <typename T, typename L, typename D, typename M>
    Kokkos::View<T*,L,D,M>
    create_view(const std::string& label, size_t size) {
      return Kokkos::View<T*,L,D,M>(label, size);
    }

#define COMPAT_INSTANT(T,D) \
    template Kokkos::View<T*,D> getKokkosViewDeepCopy<D,T>(const Teuchos::ArrayView<T>& a); \
    template Kokkos::View<const T*,D> getKokkosViewDeepCopy<D,T>(const Teuchos::ArrayView<const T>& a); \
    template void realloc<T,D>(Kokkos::View<T*,D>& v, const D::size_type size); \
    template Kokkos::View<T*,Kokkos::LayoutLeft,D,void> create_view<T,Kokkos::LayoutLeft,D,void>(const std::string& label, size_t size); \
    template Kokkos::View<T*,Kokkos::LayoutRight,D,void> create_view<T,Kokkos::LayoutRight,D,void>(const std::string& label, size_t size); \
    template Kokkos::View<T*,D,void,void> create_view<T,D,void,void>(const std::string& label, size_t size);

  } // namespace Compat
} // namespace Kokkos

#endif // KOKKOSCOMPAT_VIEW_DEF_HPP
