// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_MESH_BASE_NGPFIELDBLAS_HPP
#define STK_MESH_BASE_NGPFIELDBLAS_HPP

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/baseImpl/NgpFieldBLASImpl.hpp>

#include <complex>
#include <string>
#include <iostream>
#include <algorithm>

namespace stk {
namespace mesh {

template<class Scalar, typename EXEC_SPACE>
inline
void field_fill(const Scalar alpha, const FieldBase& field, const EXEC_SPACE& execSpace,
                bool IsDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  constexpr bool isActuallyDeviceExecSpace = !Kokkos::SpaceAccessibility<EXEC_SPACE, stk::ngp::HostExecSpace::memory_space>::accessible;
#ifdef STK_USE_DEVICE_MESH
  constexpr bool operateOnDevice = isActuallyDeviceExecSpace;
#else
  constexpr bool operateOnDevice = false;
#endif

  field.clear_sync_state();

  if constexpr (operateOnDevice) {
    NgpField<Scalar>& ngpField = get_updated_ngp_field<Scalar>(field);
    impl::field_fill_no_sync_or_mark(alpha, ngpField, execSpace);
  }
  else {
    stk::mesh::field_fill(alpha, field);
  }

  const bool markModifiedOnDevice = isActuallyDeviceExecSpace || IsDeviceExecSpaceUserOverride;
  if (markModifiedOnDevice) {
    field.modify_on_device();
  }
  else {
    field.modify_on_host();
  }
}

template<class Scalar, typename EXEC_SPACE>
inline
void field_fill(const Scalar alpha, const FieldBase& field, const Selector& selector, const EXEC_SPACE& execSpace,
                bool IsDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  constexpr bool isActuallyDeviceExecSpace = !Kokkos::SpaceAccessibility<EXEC_SPACE, stk::ngp::HostExecSpace::memory_space>::accessible;
#ifdef STK_USE_DEVICE_MESH
  constexpr bool operateOnDevice = isActuallyDeviceExecSpace;
#else
  constexpr bool operateOnDevice = false;
#endif

  field.clear_sync_state();

  if constexpr (operateOnDevice) {
    NgpField<Scalar>& ngpField = get_updated_ngp_field<Scalar>(field);
    impl::field_fill_no_sync_or_mark(alpha, ngpField, selector, execSpace);
  }
  else {
    stk::mesh::field_fill(alpha, field, selector);
  }

  const bool markModifiedOnDevice = isActuallyDeviceExecSpace || IsDeviceExecSpaceUserOverride;
  if (markModifiedOnDevice) {
    field.modify_on_device();
  }
  else {
    field.modify_on_host();
  }
}

template<typename EXEC_SPACE>
inline
void field_copy(const FieldBase& xField, const FieldBase& yField, const EXEC_SPACE& execSpace,
                bool IsDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  constexpr bool isActuallyDeviceExecSpace = !Kokkos::SpaceAccessibility<EXEC_SPACE, stk::ngp::HostExecSpace::memory_space>::accessible;
#ifdef STK_USE_DEVICE_MESH
  constexpr bool operateOnDevice = isActuallyDeviceExecSpace;
#else
  constexpr bool operateOnDevice = false;
#endif

  yField.clear_sync_state();

  if constexpr (operateOnDevice) {
    xField.sync_to_device();
    impl::field_copy_no_sync_or_mark(xField, yField, execSpace);
  }
  else {
    xField.sync_to_host();
    stk::mesh::field_copy(xField, yField);
  }

  yField.clear_sync_state();
  const bool markModifiedOnDevice = isActuallyDeviceExecSpace || IsDeviceExecSpaceUserOverride;
  if (markModifiedOnDevice) {
    yField.modify_on_device();
  }
  else {
    yField.modify_on_host();
  }
}

template<typename EXEC_SPACE>
inline
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector, const EXEC_SPACE& execSpace,
                bool IsDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  constexpr bool isActuallyDeviceExecSpace = !Kokkos::SpaceAccessibility<EXEC_SPACE, stk::ngp::HostExecSpace::memory_space>::accessible;
#ifdef STK_USE_DEVICE_MESH
  constexpr bool operateOnDevice = isActuallyDeviceExecSpace;
#else
  constexpr bool operateOnDevice = false;
#endif

  yField.clear_sync_state();

  if constexpr (operateOnDevice) {
    xField.sync_to_device();
    impl::field_copy_no_sync_or_mark(xField, yField, selector, execSpace);
  }
  else {
    xField.sync_to_host();
    stk::mesh::field_copy(xField, yField, selector);
  }

  yField.clear_sync_state();
  const bool markModifiedOnDevice = isActuallyDeviceExecSpace || IsDeviceExecSpaceUserOverride;
  if (markModifiedOnDevice) {
    yField.modify_on_device();
  }
  else {
    yField.modify_on_host();
  }
}

} // mesh
} // stk

#endif // STK_MESH_BASE_NGPFIELDBLAS_HPP

