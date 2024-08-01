// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_GETENTRYONHOST_HPP
#define TPETRA_DETAILS_GETENTRYONHOST_HPP

/// \file Tpetra_Details_getEntryOnHost.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::getEntryOnHost.
/// \warning The contents of this file are implementation details of
///   Tpetra.  We make no promises of backwards compatibility.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

template<class ViewType,
         class IndexType>
typename ViewType::non_const_value_type
getEntryOnHost (const ViewType& x,
                const IndexType ind)
{
  using execution_space = typename ViewType::execution_space;
  static_assert (ViewType::rank == 1, "x must be a rank-1 Kokkos::View.");
  // Get a 0-D subview of the entry of the array, and copy to host scalar.
  typename ViewType::non_const_value_type val;
  // DEEP_COPY REVIEW - DEVICE-TO-VALUE
  Kokkos::deep_copy(execution_space(), val, Kokkos::subview(x, ind));
  return val;
}

template<class ViewType,
         class IndexType>
typename ViewType::HostMirror::const_type
getEntriesOnHost (const ViewType& x,
                  const IndexType ind,
                  const int count)
{
  static_assert (ViewType::rank == 1, "x must be a rank-1 Kokkos::View.");
  // Get a 1-D subview of the entry of the array, and copy to host.
  auto subview = Kokkos::subview(x, Kokkos::make_pair(ind, ind + count));
  return Kokkos::create_mirror_view_and_copy(typename ViewType::HostMirror::memory_space(), subview);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GETENTRYONHOST_HPP
