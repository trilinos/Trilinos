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

#ifndef KOKKOSCOMPAT_VIEW_HPP
#define KOKKOSCOMPAT_VIEW_HPP

/// \file KokkosCompat_View.hpp
/// \brief Compatibility between Kokkos Array View and Teuchos memory management classes.
///
/// \warning This file and everything in it should be considered an
///   implementation detail of Tpetra.  We make no promises of
///   backwards compatibility for any interface in this file, nor do
///   we even promise that this header file will continue to exist.

#include <TeuchosKokkosCompat_config.h>
#include <Kokkos_Core.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Kokkos {
  namespace Compat {
    // Convert Kokkos::View to Teuchos::ArrayView
    template <typename ViewType>
    Teuchos::ArrayView<typename ViewType::value_type>
    getArrayView (const ViewType& a) {
      const auto size = a.size ();
      return Teuchos::ArrayView<typename ViewType::value_type>
        (size == 0 ? nullptr : a.data (), size);
    }
    template <typename ViewType>
    Teuchos::ArrayView<const typename ViewType::value_type>
    getConstArrayView (const ViewType& a) {
      const auto size = a.size ();
      return Teuchos::ArrayView<const typename ViewType::value_type>
        (size == 0 ? nullptr : a.data (), size);
    }

   // Convert Teuchos::ArrayView to Kokkos::View through deep_copy
    template <typename D, typename T>
    Kokkos::View<T*,D>
    getKokkosViewDeepCopy (const Teuchos::ArrayView<T>& a)
    {
      typedef typename std::conditional<
        SpaceAccessibility< D, Kokkos::HostSpace>::accessible,
        typename D::execution_space, Kokkos::HostSpace>::type
        HostDevice;
      typedef Kokkos::View<T*, D> view_type;
      typedef Kokkos::View<T*, typename view_type::array_layout, HostDevice,
                           Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
      if (a.size () == 0) {
        return view_type ();
      }
      view_type v ("", a.size ());
      unmanaged_host_view_type hv (a.getRawPtr (), a.size ());
      Kokkos::deep_copy (v, hv);
      return v;
    }

    template <typename D, typename T>
    Kokkos::View<const T*,D>
    getKokkosViewDeepCopy(const Teuchos::ArrayView<const T>& a)
    {
      typedef typename std::conditional<
        SpaceAccessibility< D, Kokkos::HostSpace>::accessible,
        typename D::execution_space, Kokkos::HostSpace>::type
        HostDevice;
      typedef Kokkos::View<T*, D>  view_type;
      typedef Kokkos::View<const T*, typename view_type::array_layout, HostDevice,
                           Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
      if (a.size () == 0) {
        return view_type ();
      }
      view_type v ("", a.size ());
      unmanaged_host_view_type hv (a.getRawPtr (), a.size ());
      Kokkos::deep_copy (v, hv);
      return v;
    }

    // Reallocate Kokkos::View to a new size
    // Currently this is only here to hide the Cuda device
    template <typename T, typename D>
    void
    realloc(Kokkos::View<T*,D>& v, const typename D::size_type size) {
      Kokkos::realloc(v,size);
    }

    template <class ViewType>
    ViewType
    create_view (const std::string& label, size_t size) {
      static_assert(Kokkos::is_view<ViewType>::value==true,"Kokkos::Compat::create_view() called with non-view argument.");
      static_assert(ViewType::rank==1,"Kokkos::Compat::create_view() called with non-rank-1 view argument.");
      return ViewType (label, size);
    }

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
      return Teuchos::arcp (view.span () == 0 ? nullptr : view.data (),
                            0, view.span (),
                            deallocator(view), false);
    }

    // Create a "persisting view" from a Kokkos::View
    template <typename ViewType>
    Teuchos::ArrayRCP<typename ViewType::value_type>
    persistingView(
      const ViewType& view,
      typename Teuchos::ArrayRCP<typename ViewType::value_type>::size_type offset,
      typename Teuchos::ArrayRCP<typename ViewType::value_type>::size_type size)
    {
      // This check is related to #4234.
      typename ViewType::value_type* const ptr =
        (size == 0 && offset == 0) ? nullptr : (view.data () + offset);
      return Teuchos::arcp (ptr, 0, size, deallocator (view), false);
    }

    template <class ViewType, typename Ordinal>
    ViewType
    subview_range (const ViewType& view,
                   const Ordinal begin,
                   const Ordinal end)
    {
      static_assert(Kokkos::is_view<ViewType>::value==true,"Kokkos::Compat::subview_range() called with non-view argument.");
      static_assert(ViewType::rank==1,"Kokkos::Compat::subview_range() called with non-rank-1 view argument.");
      return Kokkos::subview (view, std::make_pair (begin, end));
    }

    template <class ViewType, typename Ordinal>
    ViewType
    subview_offset (const ViewType& view,
                    const Ordinal offset,
                    const Ordinal size)
    {
      static_assert(Kokkos::is_view<ViewType>::value==true,"Kokkos::Compat::subview_offset() called with non-view argument.");
      static_assert(ViewType::rank==1,"Kokkos::Compat::subview_offset() called with non-rank-1 view argument.");
      return Kokkos::subview (view, std::make_pair (offset, offset+size));
    }

    template <class DstViewType, class SrcViewType,
              typename Ordinal>
    void
    deep_copy_range (const DstViewType& dst,
                     const SrcViewType& src,
                     const Ordinal dst_begin,
                     const Ordinal src_begin,
                     const Ordinal src_end)
    {
      static_assert(Kokkos::is_view<DstViewType>::value==true,"Kokkos::Compat::deep_copy_range() called with non-view argument.");
      static_assert(Kokkos::is_view<SrcViewType>::value==true,"Kokkos::Compat::deep_copy_range() called with non-view argument.");
      static_assert(DstViewType::rank==1 && SrcViewType::rank==1,"Kokkos::Compat::deep_copy_range() called with non-rank-1 view argument.");
      const Ordinal size = src_end - src_begin;
      const Ordinal dst_end = dst_begin + size;
      DstViewType dst_sub = Kokkos::subview(
        dst, std::make_pair (dst_begin, dst_end));
      SrcViewType src_sub = Kokkos::subview(
        src, std::make_pair (src_begin, src_end));
      Kokkos::deep_copy(dst_sub, src_sub);
    }

    template <class DstViewType, class SrcViewType,
              typename Ordinal>
    void
    deep_copy_offset (const DstViewType& dst,
                      const SrcViewType& src,
                      const Ordinal dst_offset,
                      const Ordinal src_offset,
                      const Ordinal size)
    {
      static_assert(Kokkos::is_view<DstViewType>::value==true,"Kokkos::Compat::deep_copy_offset() called with non-view argument.");
      static_assert(Kokkos::is_view<SrcViewType>::value==true,"Kokkos::Compat::deep_copy_offset() called with non-view argument.");
      static_assert(DstViewType::rank==1 && SrcViewType::rank==1,"Kokkos::Compat::deep_copy_offset() called with non-rank-1 view argument.");
      const Ordinal dst_end = dst_offset + size;
      const Ordinal src_end = src_offset + size;
      DstViewType dst_sub = Kokkos::subview(
        dst, std::make_pair (dst_offset, dst_end));
      SrcViewType src_sub = Kokkos::subview(
        src, std::make_pair (src_offset, src_end));
      Kokkos::deep_copy(typename SrcViewType::execution_space(), dst_sub, src_sub);
    }

    template <class ViewType>
    typename ViewType::const_type
    create_const_view(const ViewType& view) {
      static_assert(Kokkos::is_view<ViewType>::value==true,"Kokkos::Compat::create_const_view() called with non-view argument.");
      return view;
    }

  } // namespace Compat
} // namespace Kokkos

#endif // KOKKOSCOMPAT_VIEW_HPP
