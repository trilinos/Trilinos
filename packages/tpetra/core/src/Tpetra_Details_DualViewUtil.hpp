// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_DUALVIEWUTIL_HPP
#define TPETRA_DETAILS_DUALVIEWUTIL_HPP

#include "TpetraCore_config.h"
#include "Kokkos_DualView.hpp"
#include "Teuchos_ArrayView.hpp"
#include <ostream>
#include <string>

//! Namespace for Tpetra classes and methods
namespace Tpetra {

/// \brief Namespace for Tpetra implementation details.
/// \warning Do NOT rely on the contents of this namespace.
namespace Details {

/// \brief Use in place of the string label as the first argument of
///   Kokkos::View's constructor, in case you want to allocate without
///   initializing.
auto view_alloc_no_init (const std::string& label) ->
  decltype (Kokkos::view_alloc (label, Kokkos::WithoutInitializing));

/// \brief Initialize \c dv such that its host View is \c hostView.
///
/// This shallow copies the host View into the output DualView,
/// and syncs the output DualView to device.
template<class ElementType, class DeviceType>
void
makeDualViewFromOwningHostView
  (Kokkos::DualView<ElementType*, DeviceType>& dv,
   const typename Kokkos::DualView<ElementType*, DeviceType>::t_host& hostView)
{
  using execution_space = typename DeviceType::execution_space;
  using dual_view_type = Kokkos::DualView<ElementType*, DeviceType>;

  if (dv.extent (0) == hostView.extent (0)) {
    // We don't need to reallocate the device View.
    dv.clear_sync_state ();
    dv.modify_host ();
    dv.h_view = hostView;
    dv.sync_device ();
  }
  else {
    auto devView = Kokkos::create_mirror_view (DeviceType (), hostView);
    // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
    Kokkos::deep_copy (execution_space(), devView, hostView);
    dv = dual_view_type (devView, hostView);
  }
}

template<class ElementType, class DeviceType>
void
makeDualViewFromArrayView (Kokkos::DualView<ElementType*, DeviceType>& dv,
                           const Teuchos::ArrayView<const ElementType>& av,
                           const std::string& label)
{
  using execution_space = typename DeviceType::execution_space;
  using dual_view_type = Kokkos::DualView<ElementType*, DeviceType>;
  using host_view_type = typename dual_view_type::t_host;
  using const_host_view_type = typename host_view_type::const_type;

  const auto size = av.size ();
  const ElementType* ptr = (size == 0) ? nullptr : av.getRawPtr ();
  const_host_view_type inView (ptr, size);
  host_view_type hostView (view_alloc_no_init (label), size);
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  Kokkos::deep_copy (execution_space(), hostView, inView);

  makeDualViewFromOwningHostView (dv, hostView);
}

template<class ElementType, class DeviceType>
void
makeDualViewFromVector (Kokkos::DualView<ElementType*, DeviceType>& dv,
                        const std::vector<ElementType>& vec,
                        const std::string& label)
{
  using dual_view_type = Kokkos::DualView<ElementType*, DeviceType>;
  using execution_space = typename DeviceType::execution_space;
  using host_view_type = typename dual_view_type::t_host;
  using const_host_view_type = typename host_view_type::const_type;

  const auto size = vec.size ();
  const ElementType* ptr = (size == 0) ? nullptr : vec.data ();
  const_host_view_type inView (ptr, size);
  host_view_type hostView (view_alloc_no_init (label), size);
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  Kokkos::deep_copy (execution_space(), hostView, inView);

  makeDualViewFromOwningHostView (dv, hostView);
}

template<class ElementType, class DeviceType>
void
printDualView (std::ostream& out,
               const Kokkos::DualView<ElementType*, DeviceType>& dv,
               const std::string& name)
{
  out << name << ": ";
  const size_t size = size_t (dv.extent (0));
  const auto hostView = dv.view_host ();

  out << "[";
  for (size_t k = 0; k < size; ++k) {
    out << hostView[k];
    if (k + size_t (1) < size) {
      out << ",";
    }
  }
  out << "]";
}

} // namespace Details

} // namespace Tpetra

#endif // TPETRA_DETAILS_DUALVIEWUTIL_HPP
