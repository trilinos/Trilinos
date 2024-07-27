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
/// \param dv [in/out] The DualView to reallocate, if needed.
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
reallocDualViewIfNeeded (Kokkos::DualView<ValueType*, DeviceType>& dv,
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
  if (curSize == newSize) {
    return false; // did not reallocate
  }
  else if (curSize < newSize) { // too small; need to reallocate
    if (needFenceBeforeRealloc) {
      execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
    }
    dv = dual_view_type (); // free first, in order to save memory
    // If current size is 0, the DualView's Views likely lack a label.
    dv = dual_view_type (curSize == 0 ? newLabel : dv.d_view.label (), newSize);
    return true; // we did reallocate
  }
  else {
    if (newSize == 0) { // special case: realloc to 0 means always do it
      if (needFenceBeforeRealloc) {
        execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
      }
      // If current size is 0, the DualView's Views likely lack a label.
      dv = dual_view_type (curSize == 0 ? newLabel : dv.d_view.label (), 0);
      return true; // we did reallocate
    }
    // Instead of writing curSize >= tooBigFactor * newSize, express
    // via division to avoid overflow (for very large right-hand side).
    // We've already tested whether newSize == 0, so this is safe.
    else if (curSize / newSize >= tooBigFactor) {
      // The allocation is much too big, so free it and reallocate
      // to the new, smaller size.
      if (needFenceBeforeRealloc) {
        execution_space().fence (); // keep this fence to respect needFenceBeforeRealloc
      }
      dv = dual_view_type (); // free first, in order to save memory
      // If current size is 0, the DualView's Views likely lack a label.
      dv = dual_view_type (curSize == 0 ? newLabel : dv.d_view.label (), newSize);
      return true; // we did reallocate
    }
    else {
      dv.d_view = Kokkos::subview (dv.d_view, range_type (0, newSize));
      dv.h_view = Kokkos::subview (dv.h_view, range_type (0, newSize));
      return false; // we did not reallocate
    }
  }
}

/// \brief Like above, but with <tt>std::string</tt> label argument.
template<class ValueType, class DeviceType>
bool
reallocDualViewIfNeeded (Kokkos::DualView<ValueType*, DeviceType>& exports,
                         const size_t newSize,
                         const std::string& newLabel,
                         const size_t tooBigFactor = 2,
                         const bool needFenceBeforeRealloc = true)
{
  return reallocDualViewIfNeeded<ValueType, DeviceType> (exports, newSize,
                                                         newLabel.c_str (),
                                                         tooBigFactor,
                                                         needFenceBeforeRealloc);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_REALLOCDUALVIEWIFNEEDED_HPP
