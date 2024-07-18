// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_GET1DCONSTVIEW_HPP
#define TPETRA_DETAILS_GET1DCONSTVIEW_HPP

///
/// \file Tpetra_Details_get1DConstView.hpp
/// \brief Create a Kokkos::View from a raw host array.
///

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Util.hpp"
#include "Kokkos_DualView.hpp"
#include <Teuchos_Array.hpp>
#include <utility>

namespace Tpetra {
namespace Details {

// mfh 28 Apr 2016: Sometimes we have a raw host array, and we need
// to make a Kokkos::View out of it that lives in a certain memory
// space.  We don't want to make a deep copy of the input array if
// we don't need to, but if the memory spaces are different, we need
// to.  The following code does that.  The struct is an
// implementation detail, and the "free" function
// get1DConstViewOfUnmanagedArray is the interface to call.

template<class ST, class DT,
         const bool outputIsHostMemory =
           std::is_same<typename DT::memory_space, Kokkos::HostSpace>::value>
struct Get1DConstViewOfUnmanagedHostArray {};

template<class ST, class DT>
struct Get1DConstViewOfUnmanagedHostArray<ST, DT, true> {
  typedef Kokkos::View<const ST*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> output_view_type;

  static output_view_type
    getView (const char /* label */ [], const ST* x_raw, const size_t x_len)
    {
      // We can return the input array, wrapped as an unmanaged View.
      // Ignore the label, since unmanaged Views don't have labels.
      return output_view_type (x_raw, x_len);
    }
};

template<class ST, class DT>
struct Get1DConstViewOfUnmanagedHostArray<ST, DT, false> {
  typedef Kokkos::View<const ST*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> input_view_type;
  typedef Kokkos::View<const ST*, DT> output_view_type;

  static output_view_type
    getView (const char label[], const ST* x_raw, const size_t x_len)
    {
      input_view_type x_in (x_raw, x_len);
      // The memory spaces are different, so we have to create a new
      // View which is a deep copy of the input array.
      //
      // FIXME (mfh 28 Apr 2016) This needs to be converted to
      // std::string, else the compiler can't figure out what
      // constructor we're calling.
      Kokkos::View<ST*, DT> x_out (std::string (label), x_len);
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (x_out, x_in);
      return x_out;
    }
};

template<class ST, class DT>
  typename Get1DConstViewOfUnmanagedHostArray<ST, DT>::output_view_type
get1DConstViewOfUnmanagedHostArray (const char label[], const ST* x_raw, const size_t x_len)
{
  return Get1DConstViewOfUnmanagedHostArray<ST, DT>::getView (label, x_raw, x_len);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GET1DCONSTVIEW_HPP
