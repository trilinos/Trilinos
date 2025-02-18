// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_CASTAWAYCONSTDUALVIEW_HPP
#define TPETRA_DETAILS_CASTAWAYCONSTDUALVIEW_HPP

/// \file Tpetra_Details_castAwayConstDualView.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::castAwayConstDualView, an implementation detail
///   of Tpetra.
///
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  They may disappear or change at any time.

#include "Kokkos_DualView.hpp"

namespace Tpetra {
namespace Details {

/// \brief Cast away const-ness of a 1-D Kokkos::DualView.
///
/// Kokkos::DualView<const ValueType*, DeviceType> forbids sync, at
/// run time.  If we want to sync it, we have to cast away const.
template<class ValueType, class DeviceType>
Kokkos::DualView<ValueType*, DeviceType>
castAwayConstDualView (const Kokkos::DualView<const ValueType*, DeviceType>& input_dv)
{
  typedef Kokkos::DualView<const ValueType*, DeviceType> input_dual_view_type;
  typedef typename input_dual_view_type::t_dev::non_const_type out_dev_view_type;
  typedef typename input_dual_view_type::t_host::non_const_type out_host_view_type;

  out_dev_view_type output_view_dev
    (const_cast<ValueType*> (input_dv.d_view.data ()),
     input_dv.d_view.extent (0));
  out_host_view_type output_view_host
    (const_cast<ValueType*> (input_dv.h_view.data ()),
     input_dv.h_view.extent (0));

  Kokkos::DualView<ValueType*, DeviceType> output_dv(output_view_dev,output_view_host);
  if(input_dv.need_sync_host()) output_dv.modify_device();
  if(input_dv.need_sync_device()) output_dv.modify_host();
  return output_dv;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CASTAWAYCONSTDUALVIEW_HPP
