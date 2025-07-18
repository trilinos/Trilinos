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
/// \param parent_view [in/out] Parent view of the DualView current_view.
/// \param current_view [in/out] The DualView to reallocate, if needed.  This
///   can be set to the same thing as parent_view, but then you lose the
///   memory advantages of specifying the parent_view.
/// \param newSize [in] New (requested) size of the DualView.
/// \param newLabel [in] New label for the DualView; only used if
///   reallocating.
/// \param needFenceBeforeRealloc [in] Whether we need to execute a
///   fence before reallocation (see below).  The fence will only
///   happen if this function needs to reallocate.
///
/// If \c current_view is too small, reallocate it to the requested size.  If it
/// is too large, and at least tooBigFactor times bigger than it needs
/// to be, free it and reallocate to the size we need, if -DTpetra_ENABLE_Legacy_Import_Shrink=ON
/// is set in the cmake configuration to save space.  Otherwise, just set it to a subview of itself, so
/// that the size is correct.
///
/// \return Whether we actually reallocated.  If we did reallocate,
///   the function promises to fence before returning.  "Fence" means
///   <tt>DeviceType::execution_space().fence()</tt>.
template<class ValueType, class DeviceType>
bool
reallocDualViewIfNeeded (Kokkos::DualView<ValueType*, DeviceType>& parent_view,
                         Kokkos::DualView<ValueType*, DeviceType>& current_view,
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

  const size_t curSize = static_cast<size_t> (current_view.extent (0));
  const size_t parentSize = static_cast<size_t> (parent_view.extent (0));

  if (curSize == newSize) {
    // The right size, do not reallocate
    return false;
  }
  else if (curSize < newSize && newSize <= parentSize) {
    // curr is too small, but parent is big enough
    // current_view is too small, but parent is big enough
    current_view = dual_view_type ();
    auto d_view = Kokkos::subview (parent_view.view_device(), range_type (0, newSize));
    auto h_view = Kokkos::subview (parent_view.view_host(), range_type (0, newSize));
    current_view = Kokkos::DualView<ValueType*, DeviceType>(d_view, h_view);
    return false; // we did not reallocate
  }
  else if (curSize < newSize && parentSize < newSize) {
    // curr and parent are both too small, so realloc
    if (needFenceBeforeRealloc) {
      execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
    }
    current_view = dual_view_type (); // free first, in order to save memory

    // If current size is 0, the DualView's Views likely lack a label.
    parent_view = dual_view_type (curSize == 0 ? newLabel : current_view.view_device().label (), newSize);
    current_view = parent_view;
    return true; // we did reallocate
  }
  else if (newSize == 0) {
    // New size is zero
    // We've asked for a size zero vector, so we make current_view a zero-length view
    if (needFenceBeforeRealloc) {
      execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
    }
    // If current size is 0, the DualView's Views likely lack a label.
    current_view = dual_view_type (newLabel, 0);
    return true;
  }
#ifdef HAVE_TPETRA_LEGACY_IMPORT_SHRINK
  else if (curSize / newSize >= tooBigFactor) {
    // New size is too big so we shrink (and shrinking is enabled)
    // current_view is bigger than we need, so check tooBigFactor to see if we shrink

    // Instead of writing curSize >= tooBigFactor * newSize, express
    // via division to avoid overflow (for very large right-hand side).
    // We've already tested whether newSize == 0, so this is safe.

    // The allocation is much too big, so free it and reallocate
    // to the new, smaller size.
    if (needFenceBeforeRealloc) {
      execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
    }
    current_view = dual_view_type (); // free first, in order to save memory
    parent_view = dual_view_type(); // free first, in order to save memory
    // If current size is 0, the DualView's Views likely lack a label.
    parent_view = dual_view_type (curSize == 0 ? newLabel : current_view.view_device().label (), newSize);
    current_view = parent_view;
    return true; // we did reallocate
  }
#endif
  else {
    // Either we're not shrinking ever or we're not shrinking because the array isn't too big.
    auto d_view = Kokkos::subview (parent_view.view_device(), range_type (0, newSize));
    auto h_view = Kokkos::subview (parent_view.view_host(), range_type (0, newSize));
    current_view = Kokkos::DualView<ValueType*, DeviceType>(d_view, h_view);
    return false; // we did not reallocate
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
