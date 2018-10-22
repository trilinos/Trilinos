// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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

  Kokkos::DualView<ValueType*, DeviceType> output_dv;
  output_dv.d_view = output_view_dev;
  output_dv.h_view = output_view_host;
  output_dv.modified_device = input_dv.modified_device;
  output_dv.modified_host = input_dv.modified_host;
  return output_dv;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CASTAWAYCONSTDUALVIEW_HPP
