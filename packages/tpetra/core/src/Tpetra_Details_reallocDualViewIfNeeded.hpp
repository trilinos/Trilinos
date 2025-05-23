// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_REALLOCDUALVIEWIFNEEDED_HPP
#define TPETRA_DETAILS_REALLOCDUALVIEWIFNEEDED_HPP

/// \file Tpetra_Details_reallocDualViewIfNeeded.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::reallocDualViewIfNeeded, an implementation
///   detail of Tpetra.
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra.  DO NOT DEPEND ON THEM.

#include "Tpetra_Details_Profiling.hpp"
#include "Kokkos_DualView.hpp"

namespace Tpetra {
namespace Details {

/// \brief Reallocate the DualView in/out argument, if needed.
///
/// \param parent_view [in/out] Parent view of the DualView dv.
/// \param dv [in/out] The DualView to reallocate, if needed.  This
///   can be set to the same thing as parent_view, but then you lose the
///   memory advantages of specifying the parentview.
/// \param newSize [in] New (requested) size of the DualView.
/// \param newLabel [in] New label for the DualView; only used if
///   reallocating.
/// \param tooBigFactor [in] Factor for deciding whether to free and
///   reallocate, or just take a subview, if \c dv is too big.  Taking
///   a subview avoids reallocation, which is expensive for some
///   memory spaces.
/// \param needFenceBeforeRealloc [in] Whether we need to execute a
///   fence before reallocation (see below).  The fence will only
///   happen if this function needs to reallocate.
///
/// If \c dv is too small, reallocate it to the requested size.  If it
/// is too large, and at least tooBigFactor times bigger than it needs
/// to be, free it and reallocate to the size we need, in order to
/// save space.  Otherwise, just set it to a subview of itself, so
/// that the size is correct.
///
/// \return Whether we actually reallocated.  If we did reallocate,
///   the function promises to fence before returning.  "Fence" means
///   <tt>DeviceType::execution_space().fence()</tt>.
template<class ValueType, class DeviceType>
bool
reallocDualViewIfNeeded (Kokkos::DualView<ValueType*, DeviceType>& parent_view,
                         Kokkos::DualView<ValueType*, DeviceType>& dv,
                         const size_t newSize,
                         const char newLabel[],
                         const size_t tooBigFactor = 2,
                         const bool needFenceBeforeRealloc = true)
{
  typedef typename DeviceType::execution_space execution_space;
  typedef Kokkos::DualView<ValueType*, DeviceType> dual_view_type;
  typedef Kokkos::pair<size_t, size_t> range_type;

  // Profiling this matters, because GPU allocations can be expensive.
  using Tpetra::Details::ProfilingRegion;
  ProfilingRegion region ("Tpetra::Details::reallocDualViewIfNeeded");

  const size_t curSize = static_cast<size_t> (dv.extent (0));
  const size_t parentSize = static_cast<size_t> (parent_view.extent (0));

  if (curSize == newSize) {
    // The right size, do not reallocate
    std::cout<<"reallocDualViewIfNeeded: Correct size = "<<newSize<<std::endl;
    return false;
  }
  else if (curSize < newSize && newSize <= parentSize) {
    std::cout<<"reallocDualViewIfNeeded: curr is too small, but parent is big enough "<<std::endl;
    // dv is too small, but parent is big enough
    dv = dual_view_type ();
    auto d_view = Kokkos::subview (parent_view.view_device(), range_type (0, newSize));
    auto h_view = Kokkos::subview (parent_view.view_host(), range_type (0, newSize));
    dv = Kokkos::DualView<ValueType*, DeviceType>(d_view, h_view);
    return false; // we did not reallocate
  }
  else if (curSize < newSize && parentSize < newSize) {
    std::cout<<"reallocDualViewIfNeeded: curr and parent are too small "<<std::endl;
    // both the parent view and dv are too small
    if (needFenceBeforeRealloc) {
      execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
    }
    dv = dual_view_type (); // free first, in order to save memory

    // If current size is 0, the DualView's Views likely lack a label.
    parent_view = dual_view_type (curSize == 0 ? newLabel : dv.view_device().label (), newSize);
    dv = parent_view;
    return true; // we did reallocate
  }
  else if (newSize == 0) {
    std::cout<<"reallocDualViewIfNeeded: new size is zero "<<std::endl;
    // We've asked for a size zero vector, so we make dv a zero-length view
    if (needFenceBeforeRealloc) {
      execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
    }
    // If current size is 0, the DualView's Views likely lack a label.
    dv = dual_view_type (newLabel, 0);
    return true;
  }
  else {
#ifdef ENABLE_SHRINK
    std::cout<<"reallocDualViewIfNeeded: new size is too big so we shrink newSize = "<newSize<std::endl;
    // dv is bigger than we need, so check tooBigFactor to see if we shrink

    // Instead of writing curSize >= tooBigFactor * newSize, express
    // via division to avoid overflow (for very large right-hand side).
    // We've already tested whether newSize == 0, so this is safe.
    if (curSize / newSize >= tooBigFactor) {
      // The allocation is much too big, so free it and reallocate
      // to the new, smaller size.
      if (needFenceBeforeRealloc) {
        execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
      }
      dv = dual_view_type (); // free first, in order to save memory
      parent_view = dual_view_type(); // free first, in order to save memory
      // If current size is 0, the DualView's Views likely lack a label.
      parent_view = dual_view_type (curSize == 0 ? newLabel : dv.view_device().label (), newSize);
      dv = parent_view;
      return true; // we did reallocate
    }
    else {
#else
      std::cout<<"reallocDualViewIfNeeded: new size is too big so we should shrink (but shrink is disabled) newSize = "<<newSize<<std::endl;
#endif
      auto d_view = Kokkos::subview (parent_view.view_device(), range_type (0, newSize));
      auto h_view = Kokkos::subview (parent_view.view_host(), range_type (0, newSize));
      dv = Kokkos::DualView<ValueType*, DeviceType>(d_view, h_view);
      return false; // we did not reallocate
#ifdef ENABLE_SHRINK
    }
#endif
  }
}

/// \brief Like above, but with <tt>std::string</tt> label argument.
template<class ValueType, class DeviceType>
bool
reallocDualViewIfNeeded (Kokkos::DualView<ValueType*, DeviceType>& parent_view,
                         Kokkos::DualView<ValueType*, DeviceType>& exports,
                         const size_t newSize,
                         const std::string& newLabel,
                         const size_t tooBigFactor = 2,
                         const bool needFenceBeforeRealloc = true)
{
  return reallocDualViewIfNeeded<ValueType, DeviceType> (parent_view,
                                                         exports, newSize,
                                                         newLabel.c_str (),
                                                         tooBigFactor,
                                                         needFenceBeforeRealloc);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_REALLOCDUALVIEWIFNEEDED_HPP
