// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_CREATEMIRRORVIEW_HPP
#define TPETRA_DETAILS_CREATEMIRRORVIEW_HPP

#include "TpetraCore_config.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Kokkos_Core.hpp"
#include <memory>
#include <string>

/// \file Tpetra_Details_createMirrorView.hpp
/// \brief Functions that wrap Kokkos::create_mirror_view, in order to
///   avoid deep copies when not necessary, even if the inputs are
///   const.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {

namespace Impl {

// Implementation detail of create_mirror_view_from_raw_host_array
// (see below).
template<class ValueType,
         class OutputDeviceType,
         const bool constInput = std::is_const<ValueType>::value,
         const bool sameAsHost =
           std::is_same<Kokkos::HostSpace,
             typename OutputDeviceType::memory_space>::value>
class CreateMirrorViewFromUnmanagedHostArray {
public:
  typedef Kokkos::View<ValueType*, OutputDeviceType> output_view_type;
  typedef Kokkos::View<ValueType*,
                       typename output_view_type::array_layout,
                       Kokkos::HostSpace> input_view_type;
  static output_view_type
  doIt (ValueType* inPtr,
        const size_t inSize,
        const bool copy = true,
        const char label[] = "");
};

// Implementation detail of create_mirror_view_from_raw_host_array
// (see below).
template<class ValueType,
         class OutputDeviceType,
         const bool constInput>
class CreateMirrorViewFromUnmanagedHostArray<ValueType, OutputDeviceType, constInput, true> {
public:
  typedef Kokkos::View<ValueType*, OutputDeviceType> output_view_type;
  typedef Kokkos::View<ValueType*, typename output_view_type::array_layout,
                       Kokkos::HostSpace> input_view_type;
  static output_view_type
  doIt (ValueType* inPtr,
        const size_t inSize,
        const bool /* copy */,
        const char /* label */ [] = "")
  {
    static_assert (std::is_same<typename OutputDeviceType::memory_space,
                     Kokkos::HostSpace>::value,
                   "OutputDeviceType::memory_space must be the same as "
                   "Kokkos::HostSpace in order to use this specialization.  "
                   "Please report this bug to the Tpetra developers.");
    return output_view_type (inPtr, inSize);
  }
};

// Implementation detail of create_mirror_view_from_raw_host_array
// (see below).
template<class ValueType,
         class OutputDeviceType>
class CreateMirrorViewFromUnmanagedHostArray<ValueType, OutputDeviceType, true, false> {
public:
  typedef Kokkos::View<ValueType*, OutputDeviceType> output_view_type;
  typedef Kokkos::View<ValueType*, typename output_view_type::array_layout,
                       Kokkos::HostSpace> input_view_type;
  static output_view_type
  doIt (ValueType* inPtr,
        const size_t inSize,
        const bool copy = true,
        const char label[] = "")
  {
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    static_assert (std::is_const<ValueType>::value, "ValueType must be const "
                   "in order to use this specialization.  Please report this "
                   "bug to the Tpetra developers.");
    static_assert (! std::is_same<typename OutputDeviceType::memory_space, Kokkos::HostSpace>::value,
                   "OutputDeviceType::memory_space must not be the same as "
                   "Kokkos::HostSpace in order to use this specialization.  "
                   "Please report this bug to the Tpetra developers.");
    input_view_type inView (inPtr, inSize);
    // ValueType is const, so we have to strip away const first.
    typedef typename output_view_type::non_const_type nc_output_view_type;
    nc_output_view_type outView_nc;
    if (! copy) {
      // Label needs to be a string and not a char*, if given as an
      // argument to Kokkos::view_alloc.  This is because view_alloc
      // also allows a raw pointer as its first argument.  See
      // https://github.com/kokkos/kokkos/issues/434.
      outView_nc = nc_output_view_type (view_alloc (std::string (label)), inSize);
    }
    else {
      // No need to initialize, if we're going to copy into it anyway.
      outView_nc = nc_output_view_type (view_alloc (std::string (label), WithoutInitializing), inSize);
      // DEEP_COPY REVIEW - HOST-TO-DEVICE
      using execution_space = typename nc_output_view_type::execution_space;
      Kokkos::deep_copy (execution_space(), outView_nc, inView);
    }
    return outView_nc; // this casts back to const
  }
};

// Implementation detail of create_mirror_view_from_raw_host_array
// (see below).
template<class ValueType,
         class OutputDeviceType>
class CreateMirrorViewFromUnmanagedHostArray<ValueType, OutputDeviceType, false, false> {
public:
  typedef Kokkos::View<ValueType*, OutputDeviceType> output_view_type;
  typedef Kokkos::View<ValueType*, typename output_view_type::array_layout,
                       Kokkos::HostSpace> input_view_type;
  static output_view_type
  doIt (ValueType* inPtr,
        const size_t inSize,
        const bool copy = true,
        const char label[] = "")
  {
    typedef typename OutputDeviceType::memory_space out_mem_space;
    typedef typename OutputDeviceType::execution_space out_exec_space;
    static_assert (! std::is_const<ValueType>::value, "ValueType must not be "
                   "const in order to use this specialization.  Please report "
                   "this bug to the Tpetra developers.");
    static_assert (! std::is_same<out_mem_space, Kokkos::HostSpace>::value,
                   "OutputDeviceType::memory_space must not be the same as "
                   "Kokkos::HostSpace in order to use this specialization.  "
                   "Please report this bug to the Tpetra developers.");
    input_view_type inView (inPtr, inSize);
    output_view_type outView =
      Kokkos::create_mirror_view (out_mem_space (), inView);
    if (copy) {
      // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
      Kokkos::deep_copy (out_exec_space(), outView, inView);
    }
    return outView;
  }
};

} // namespace Impl

/// \brief Variant of Kokkos::create_mirror_view that takes a raw host
///   1-d array as input.
///
/// Given a pointer to a 1-D array in host memory, and the number of
/// entries in the array, return a Kokkos::View that lives in
/// OutputDeviceType, and that is a mirror view of the input array.
/// By default, copy the host data to the output View, if necessary.
template<class ValueType, class OutputDeviceType>
typename Impl::CreateMirrorViewFromUnmanagedHostArray<ValueType, OutputDeviceType>::output_view_type
create_mirror_view_from_raw_host_array (const OutputDeviceType& /* dev */,
                                        ValueType* inPtr,
                                        const size_t inSize,
                                        const bool copy = true,
                                        const char label[] = "")
{
  typedef Impl::CreateMirrorViewFromUnmanagedHostArray<ValueType, OutputDeviceType> impl_type;
  return impl_type::doIt (inPtr, inSize, copy, label);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CREATEMIRRORVIEW_HPP
