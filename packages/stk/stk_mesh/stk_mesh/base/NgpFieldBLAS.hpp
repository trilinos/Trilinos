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
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/baseImpl/NgpFieldBLASImpl.hpp>

#include <complex>
#include <string>
#include <iostream>
#include <algorithm>

namespace stk::mesh {

namespace impl {

#ifndef STK_HIDE_DEPRECATED_CODE // delete after June 2026
template<class EXEC_SPACE>
STK_DEPRECATED
constexpr bool operate_on_ngp_mesh()
{
#ifdef STK_USE_DEVICE_MESH

#ifdef STK_ENABLE_GPU
  constexpr bool isActuallyDeviceExecSpace = !Kokkos::SpaceAccessibility<EXEC_SPACE, stk::ngp::HostExecSpace::memory_space>::accessible;
  constexpr bool operateOnNgpMesh = isActuallyDeviceExecSpace;
#else
  constexpr bool operateOnNgpMesh = true;
#endif

#else
  constexpr bool operateOnNgpMesh = false;
#endif

  return operateOnNgpMesh;
}
#endif

}

#ifndef STK_HIDE_DEPRECATED_CODE // delete after June 2026
template<class Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, field, component, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_fill(const Scalar alpha,
                const FieldBase& field,
                int component,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, field, component, selector, execSpace);
  }
  else {
    stk::mesh::field_fill<stk::ngp::HostSpace>(alpha, field, component, selector, execSpace);
  }
}

template<class Scalar, class EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, field, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_fill(const Scalar alpha,
                const FieldBase& field,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, field, execSpace);
  }
  else {
    stk::mesh::field_fill<stk::ngp::HostSpace>(alpha, field, execSpace);
  }
}

template<class Scalar, class EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, field, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_fill(const Scalar alpha,
                const FieldBase& field,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, field, selector, execSpace);
  }
  else {
    stk::mesh::field_fill<stk::ngp::HostSpace>(alpha, field, selector, execSpace);
  }
}


template<class Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, fields, component, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_fill(const Scalar alpha,
                const std::vector<const FieldBase*>& fields,
                int component,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, fields, component, selector, execSpace);
  }
  else {
    stk::mesh::field_fill<stk::ngp::HostSpace>(alpha, fields, component, selector, execSpace);
  }
}

template<class Scalar, class EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, fields, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_fill(const Scalar alpha,
                const std::vector<const FieldBase*>& fields,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, fields, execSpace);
  }
  else {
    stk::mesh::field_fill<stk::ngp::HostSpace>(alpha, fields, execSpace);
  }
}

template<class Scalar, class EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, fields, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_fill(const Scalar alpha,
                const std::vector<const FieldBase*>& fields,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_fill<stk::ngp::DeviceSpace>(alpha, fields, selector, execSpace);
  }
  else {
    stk::mesh::field_fill<stk::ngp::HostSpace>(alpha, fields, selector, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_amax<stk::ngp::DeviceSpace>(amaxOut, xField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_amax(Scalar& amaxOut,
    const FieldBase& xField,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace, EXEC_SPACE>) )
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_amax<stk::ngp::DeviceSpace>(amaxOut, xField, execSpace);
  }
  else {
    stk::mesh::field_amax<stk::ngp::HostSpace>(amaxOut, xField, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_amax<stk::ngp::DeviceSpace>(amaxOut, xField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_amax(Scalar& amaxOut,
    const FieldBase& xField,
    const Selector& selector,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace, EXEC_SPACE>) )
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_amax<stk::ngp::DeviceSpace>(amaxOut, xField, selector, execSpace);
  }
  else {
    stk::mesh::field_amax<stk::ngp::HostSpace>(amaxOut, xField, selector, execSpace);
  }
}

template<typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_copy<stk::ngp::DeviceSpace>(xField, yField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_copy(const FieldBase& xField,
                const FieldBase& yField,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_copy<stk::ngp::DeviceSpace>(xField, yField, execSpace);
  }
  else {
    stk::mesh::field_copy<stk::ngp::HostSpace>(xField, yField, execSpace);
  }
}

template<typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_copy<stk::ngp::DeviceSpace>(xField, yField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_copy(const FieldBase& xField,
                const FieldBase& yField,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_copy<stk::ngp::DeviceSpace>(xField, yField, selector, execSpace);
  }
  else {
    stk::mesh::field_copy<stk::ngp::HostSpace>(xField, yField, selector, execSpace);
  }
}

template<class DataType, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_axpby<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_axpby(const stk::mesh::BulkData& mesh,
    const DataType alpha,
    const stk::mesh::FieldBase & xField,
    const DataType beta,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_axpby<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, selector, execSpace);
  }
  else {
    stk::mesh::field_axpby<stk::ngp::HostSpace>(alpha, xField, beta, yField, selector, execSpace);
  }
}

template<class DataType, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_axpby<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_axpby(const stk::mesh::BulkData& mesh,
    const DataType alpha,
    const stk::mesh::FieldBase & xField,
    const DataType beta,
    const stk::mesh::FieldBase & yField,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_axpby<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, execSpace);
  }
  else {
    stk::mesh::field_axpby<stk::ngp::HostSpace>(alpha, xField, beta, yField, execSpace);
  }
}

template<class Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_axpbyz<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, zField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_axpbyz(const stk::mesh::BulkData& mesh,
    const Scalar alpha,
    const stk::mesh::FieldBase & xField,
    const Scalar beta,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::FieldBase & zField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_axpbyz<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, zField, selector, execSpace);
  }
  else {
    stk::mesh::field_axpbyz<stk::ngp::HostSpace>(alpha, xField, beta, yField, zField, selector, execSpace);
  }
}

template<class Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_axpbyz<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, zField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_axpbyz(const stk::mesh::BulkData& mesh,
    const Scalar alpha,
    const stk::mesh::FieldBase & xField,
    const Scalar beta,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::FieldBase & zField,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_axpbyz<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, zField, execSpace);
  }
  else {
    stk::mesh::field_axpbyz<stk::ngp::HostSpace>(alpha, xField, beta, yField, zField, execSpace);
  }
}


template<class Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_axpy<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_axpy(const stk::mesh::BulkData& mesh,
    const Scalar alpha,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_axpy<stk::ngp::DeviceSpace>(alpha, xField, yField, selector, execSpace);
  }
  else {
    stk::mesh::field_axpy<stk::ngp::HostSpace>(alpha, xField, yField, selector, execSpace);
  }
}

template<class Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_axpy<stk::ngp::DeviceSpace>(alpha, xField, beta, yField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_axpy(const stk::mesh::BulkData& mesh,
    const Scalar alpha,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_axpy<stk::ngp::DeviceSpace>(alpha, xField, yField, execSpace);
  }
  else {
    stk::mesh::field_axpy<stk::ngp::HostSpace>(alpha, xField, yField, execSpace);
  }
}

template<typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_product<stk::ngp::DeviceSpace>(xField, yField, zField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_product(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::FieldBase & zField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_product<stk::ngp::DeviceSpace>(xField, yField, zField, selector, execSpace);
  }
  else {
    stk::mesh::field_product<stk::ngp::HostSpace>(xField, yField, zField, selector, execSpace);
  }
}

template<typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_product<stk::ngp::DeviceSpace>(xField, yField, zField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_product(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::FieldBase & zField,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_product<stk::ngp::DeviceSpace>(xField, yField, zField, execSpace);
  }
  else {
    stk::mesh::field_product<stk::ngp::HostSpace>(xField, yField, zField, execSpace);
  }
}


template<typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_scale<stk::ngp::DeviceSpace>(alpha, xField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_scale(const stk::mesh::BulkData& mesh,
    const Scalar alpha,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_scale<stk::ngp::DeviceSpace>(alpha, xField, selector, execSpace);
  }
  else {
    stk::mesh::field_scale<stk::ngp::HostSpace>(alpha, xField, selector, execSpace);
  }
}

template<typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_scale<stk::ngp::DeviceSpace>(alpha, xField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_scale(const stk::mesh::BulkData& mesh,
    const Scalar alpha,
    const stk::mesh::FieldBase & xField,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_scale<stk::ngp::DeviceSpace>(alpha, xField, execSpace);
  }
  else {
    stk::mesh::field_scale<stk::ngp::HostSpace>(alpha, xField, execSpace);
  }
}

template<typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_swap<stk::ngp::DeviceSpace>(xField, yField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_swap(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_swap<stk::ngp::DeviceSpace>(xField, yField, selector, execSpace);
  }
  else {
    stk::mesh::field_swap<stk::ngp::HostSpace>(xField, yField, selector, execSpace);
  }
}

template<typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_swap<stk::ngp::DeviceSpace>(xField, yField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_swap(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_swap<stk::ngp::DeviceSpace>(xField, yField, execSpace);
  }
  else {
    stk::mesh::field_swap<stk::ngp::HostSpace>(xField, yField, execSpace);
  }
}

template<typename ReturnT, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_dot<stk::ngp::DeviceSpace>(result, xField, yField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_dot(ReturnT& result,
                      const stk::mesh::BulkData& mesh,
                      const stk::mesh::FieldBase & xField,
                      const stk::mesh::FieldBase & yField,
                      const stk::mesh::Selector & selector,
                      const EXEC_SPACE& execSpace)
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_dot<stk::ngp::DeviceSpace>(result, xField, yField, selector, execSpace);
  }
  else {
    stk::mesh::field_dot<stk::ngp::HostSpace>(result, xField, yField, selector, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_amin<stk::ngp::DeviceSpace>(result, xField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_amin(Scalar& aminOut,
    const FieldBase& xField,
    const EXEC_SPACE& execSpace)
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_amin<stk::ngp::DeviceSpace>(aminOut, xField, execSpace);
  }
  else {
    stk::mesh::field_amin<stk::ngp::HostSpace>(aminOut, xField, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_amin<stk::ngp::DeviceSpace>(result, xField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_amin(Scalar& aminOut,
    const FieldBase& xField,
    const Selector& selector,
    const EXEC_SPACE& execSpace)
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_amin<stk::ngp::DeviceSpace>(aminOut, xField, selector, execSpace);
  }
  else {
    stk::mesh::field_amin<stk::ngp::HostSpace>(aminOut, xField, selector, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_sum<stk::ngp::DeviceSpace>(result, xField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_asum(Scalar& aminOut,
    const FieldBase& xField,
    const EXEC_SPACE& execSpace )
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_asum<stk::ngp::DeviceSpace>(aminOut, xField, execSpace);
  }
  else {
    stk::mesh::field_asum<stk::ngp::HostSpace>(aminOut, xField, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_sum<stk::ngp::DeviceSpace>(result, xField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_asum(Scalar& aminOut,
    const FieldBase& xField,
    const Selector& selector,
    const EXEC_SPACE& execSpace)
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_asum<stk::ngp::DeviceSpace>(aminOut, xField, selector, execSpace);
  }
  else {
    stk::mesh::field_asum<stk::ngp::HostSpace>(aminOut, xField, selector, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_nrm2<stk::ngp::DeviceSpace>(result, xField, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_nrm2(Scalar& nrm2Out,
    const FieldBase& xField,
    const EXEC_SPACE& execSpace )
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_nrm2<stk::ngp::DeviceSpace>(nrm2Out, xField, execSpace);
  }
  else {
    stk::mesh::field_nrm2<stk::ngp::HostSpace>(nrm2Out, xField, execSpace);
  }
}

template <typename Scalar, typename EXEC_SPACE>
STK_DEPRECATED_MSG("Replace with code `stk::mesh::field_nrm2<stk::ngp::DeviceSpace>(result, xField, selector, execSpace)`.  execSpace is not required if you're just using the default.  To run on host, use `stk::ngp::HostSpace`")
void field_nrm2(Scalar& nrm2Out,
    const FieldBase& xField,
    const Selector& selector,
    const EXEC_SPACE& execSpace)
{
  if constexpr (impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::field_nrm2<stk::ngp::DeviceSpace>(nrm2Out, xField, selector, execSpace);
  }
  else {
    stk::mesh::field_nrm2<stk::ngp::HostSpace>(nrm2Out, xField, selector, execSpace);
  }
}
#endif

} // stk::mesh

#endif // STK_MESH_BASE_NGPFIELDBLAS_HPP

