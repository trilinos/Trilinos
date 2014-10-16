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

#ifndef KOKKOSCOMPAT_VIEW_HPP
#define KOKKOSCOMPAT_VIEW_HPP

/// \file KokkosCompat_View.hpp
/// \brief Compatibility between Kokkos Array View and Teuchos memory management classes.
///
/// \warning This file and everything in it should be considered an
///   implementation detail of Tpetra.  We make no promises of
///   backwards compatibility for any interface in this file, nor do
///   we even promise that this header file will continue to exist.

#include "KokkosCompat_config.h"

// KokkosCore device types
#include "Kokkos_Core.hpp"

#include "Teuchos_ArrayView.hpp"

namespace Kokkos {
  namespace Compat {

    // Convert Kokkos::View to Teuchos::ArrayView
    template <typename ViewType>
    Teuchos::ArrayView<typename ViewType::value_type>
    getArrayView(const ViewType& a) {
      return Teuchos::ArrayView<typename ViewType::value_type>(
        a.ptr_on_device(), a.size());
    }
    template <typename ViewType>
    Teuchos::ArrayView<const typename ViewType::value_type>
    getConstArrayView(const ViewType& a) {
      return Teuchos::ArrayView<const typename ViewType::value_type>(
        a.ptr_on_device(), a.size());
    }

   // Convert Teuchos::ArrayView to Kokkos::View through deep_copy
    template <typename D, typename T>
    Kokkos::View<T*,D>
    getKokkosViewDeepCopy(const Teuchos::ArrayView<T>& a);

    template <typename D, typename T>
    Kokkos::View<const T*,D>
    getKokkosViewDeepCopy(const Teuchos::ArrayView<const T>& a);

    // Reallocate Kokkos::View to a new size
    // Currently this is only here to hide the Cuda device
    template <typename T, typename D>
    void
    realloc(Kokkos::View<T*,D>& v, const typename D::size_type size);

    template <typename T, typename L, typename D, typename M>
    Kokkos::View<T*,L,D,M>
    create_view(const std::string& label, size_t size);

    // Custom deallocator for Teuchos::ArrayRCP.  It doesn't actually
    // deallocate the ArrayRCP's data (which belongs to the View passed
    // into Deallocator's constructor, not to the ArrayRCP).  Instead,
    // it just keeps a reference to the View.  That way, the data won't
    // go away until the View's reference count goes to zero.
    template<class ViewType>
    class Deallocator {
    public:
      typedef ViewType view_type;
      typedef typename ViewType::value_type ptr_t;
      typedef typename view_type::value_type value_type;

      // Constructor takes the View that owns the memory.
      Deallocator (const ViewType& view__) : view_ (view__) {}

      // "Deallocation function" doesn't actually deallocate its input
      // pointer; the View is responsible for deallocation of its
      // memory.
      void free (value_type*) {}
      view_type view() { return view_;}
    private:
      view_type view_; // View that owns the memory

    };

    // Create deallocator for a given view
    template <class ViewType>
    Deallocator<ViewType> deallocator(const ViewType& view) {
      return Deallocator<ViewType>(view);
    }

    // Create a "persisting view" from a Kokkos::View
    template <typename ViewType>
    Teuchos::ArrayRCP<typename ViewType::value_type>
    persistingView(const ViewType& view) {
      // mfh (05 Sep 2013) We use a custom deallocator that keeps the
      // view, but otherwise does nothing.  It doesn't deallocate the
      // view; it just maintains the reference count for the lifetime
      // of the deallocator object (and therefore, for the lifetime of
      // the ArrayRCP).  It seems weird that this is a nonowning
      // ArrayRCP (the last argument is false).  However, that's OK,
      // because (a) the ArrayRCP (actually, the underlying
      // RCPNodeTmpl; see Teuchos_RCPNode.hpp) keeps the deallocator
      // object regardless of whether the ArrayRCP is owning, and (b)
      // if we make this an owning ArrayRCP, the RCPNode tracking
      // system (in a debug build: Teuchos_ENABLE_DEBUG:BOOL=ON) will
      // report problems, because multiple owning ArrayRCPs have the
      // same raw pointer but not the same underlying RCPNode.  That
      // would be a problem if those ArrayRCPs shared memory that was
      // allocated (e.g.,) using new or malloc, but it's not a problem
      // here, because the custom deallocator does not free anything.
      // Nevertheless, it's better not to trouble the tracking system.
      return Teuchos::arcp(view.ptr_on_device(), 0, view.capacity(),
                           deallocator(view), false);
    }

    // Create a "persisting view" from a Kokkos::View
    template <typename ViewType>
    Teuchos::ArrayRCP<typename ViewType::value_type>
    persistingView(
      const ViewType& view,
      typename Teuchos::ArrayRCP<typename ViewType::value_type>::size_type offset,
      typename Teuchos::ArrayRCP<typename ViewType::value_type>::size_type size) {
      return Teuchos::arcp(view.ptr_on_device()+offset, 0, size,
                           deallocator(view), false);
    }

    template <typename T, typename L, typename D, typename M, typename Ordinal>
    Kokkos::View<T*,L,D,M>
    subview_range (const Kokkos::View<T*,L,D,M>& view,
                   const Ordinal begin,
                   const Ordinal end)
    {
      typedef Kokkos::View<T*,L,D,M> view_type;
      return Kokkos::subview<view_type> (view, std::make_pair (begin, end));
    }

    template <typename T, typename L, typename D, typename M, typename Ordinal>
    Kokkos::View<T*,L,D,M>
    subview_offset (const Kokkos::View<T*,L,D,M>& view,
                    const Ordinal offset,
                    const Ordinal size)
    {
      typedef Kokkos::View<T*,L,D,M> view_type;
      return Kokkos::subview<view_type> (view, std::make_pair (offset, offset+size));
    }

    template <typename DT, typename DL, typename DD, typename DM,
              typename ST, typename SL, typename SD, typename SM,
              typename Ordinal>
    void
    deep_copy_range (const Kokkos::View<DT*,DL,DD,DM>& dst,
                     const Kokkos::View<ST*,SL,SD,SM>& src,
                     const Ordinal dst_begin,
                     const Ordinal src_begin,
                     const Ordinal src_end)
    {
      typedef Kokkos::View<DT*,DL,DD,DM> dst_view_type;
      typedef Kokkos::View<ST*,SL,SD,SM> src_view_type;
      const Ordinal size = src_end - src_begin;
      const Ordinal dst_end = dst_begin + size;
      dst_view_type dst_sub = Kokkos::subview<dst_view_type>(
        dst, std::make_pair (dst_begin, dst_end));
      src_view_type src_sub = Kokkos::subview<src_view_type>(
        src, std::make_pair (src_begin, src_end));
      Kokkos::deep_copy(dst_sub, src_sub);
    }

    template <typename DT, typename DL, typename DD, typename DM,
              typename ST, typename SL, typename SD, typename SM,
              typename Ordinal>
    void
    deep_copy_offset (const Kokkos::View<DT*,DL,DD,DM>& dst,
                      const Kokkos::View<ST*,SL,SD,SM>& src,
                      const Ordinal dst_offset,
                      const Ordinal src_offset,
                      const Ordinal size)
    {
      typedef Kokkos::View<DT*,DL,DD,DM> dst_view_type;
      typedef Kokkos::View<ST*,SL,SD,SM> src_view_type;
      const Ordinal dst_end = dst_offset + size;
      const Ordinal src_end = src_offset + size;
      dst_view_type dst_sub = Kokkos::subview<dst_view_type>(
        dst, std::make_pair (dst_offset, dst_end));
      src_view_type src_sub = Kokkos::subview<src_view_type>(
        src, std::make_pair (src_offset, src_end));
      Kokkos::deep_copy(dst_sub, src_sub);
    }

    template <typename T, typename L, typename D, typename M>
    Kokkos::View<const T*, L, D, M>
    create_const_view(const Kokkos::View<T*,L,D,M>& view) {
      return view;
    }

    template <typename T, typename L, typename D, typename M>
    Kokkos::View<const T*, L, D, M>
    create_const_view(const Kokkos::View<const T*,L,D,M>& view) {
      return view;
    }

  } // namespace Compat
} // namespace Kokkos

#endif // KOKKOSCOMPAT_VIEW_HPP
