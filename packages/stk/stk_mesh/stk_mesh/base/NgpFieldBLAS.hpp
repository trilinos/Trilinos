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
#include <stk_mesh/base/FieldBase.hpp>
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
void field_fill(const Scalar alpha,
                const FieldBase& field,
                int component,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  ngp_field_blas::impl::field_fill_impl(alpha, field, component, &selector, execSpace, isDeviceExecSpaceUserOverride);
}

template<class Scalar, class EXEC_SPACE>
inline
void field_fill(const Scalar alpha,
                const FieldBase& field,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  ngp_field_blas::impl::field_fill_impl(alpha, field, -1, nullptr, execSpace, isDeviceExecSpaceUserOverride);
}

template<class Scalar, class EXEC_SPACE>
inline
void field_fill(const Scalar alpha,
                const FieldBase& field,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  ngp_field_blas::impl::field_fill_impl(alpha, field, -1, &selector, execSpace, isDeviceExecSpaceUserOverride);
}

template<typename EXEC_SPACE>
inline
void field_copy(const FieldBase& xField,
                const FieldBase& yField,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  ngp_field_blas::impl::field_copy_impl(xField, yField, nullptr, execSpace, isDeviceExecSpaceUserOverride);
}

template<typename EXEC_SPACE>
inline
void field_copy(const FieldBase& xField,
                const FieldBase& yField,
                const Selector& selector,
                const EXEC_SPACE& execSpace,
                bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  ngp_field_blas::impl::field_copy_impl(xField, yField, &selector, execSpace, isDeviceExecSpaceUserOverride);
}

template<class DataType, typename EXEC_SPACE>
inline void field_axpby(const stk::mesh::BulkData& mesh,
    const DataType alpha,
    const stk::mesh::FieldBase & xField,
    const DataType beta,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool IsDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  // y = a*x + b*y

   if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
     ngp_field_blas::impl::apply_functor_on_field<double, stk::ngp_field_blas::impl::FieldAXPBYFunctor>(
      mesh, yField, xField, yField, alpha, beta, selector);
   }
   else {
    xField.sync_to_host();
    yField.sync_to_host();
    stk::mesh::field_axpby(alpha, xField, beta, yField, selector);
  }

  yField.clear_sync_state();
  if (ngp_field_blas::impl::mark_modified_on_device(execSpace, IsDeviceExecSpaceUserOverride)) {
    yField.modify_on_device();
  }
  else {
    yField.modify_on_host();
  }
}

template<class DataType, typename EXEC_SPACE>
inline void field_axpbyz(const stk::mesh::BulkData& mesh,
    const DataType alpha,
    const stk::mesh::FieldBase & xField,
    const DataType beta,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::FieldBase & zField,
    const stk::mesh::Selector & selector,
    const EXEC_SPACE& execSpace,
    bool IsDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  // z = a*x + b*y

   if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
     ngp_field_blas::impl::apply_functor_on_field<double, stk::ngp_field_blas::impl::FieldAXPBYFunctor>(
      mesh, zField, xField, yField, alpha, beta, selector);
   }
   else {
    xField.sync_to_host();
    yField.sync_to_host();
    stk::mesh::field_copy(yField, zField, selector);
    stk::mesh::field_axpby(alpha, xField, beta, zField, selector);
  }

  zField.clear_sync_state();
  if (ngp_field_blas::impl::mark_modified_on_device(execSpace, IsDeviceExecSpaceUserOverride)) {
    zField.modify_on_device();
  }
  else {
    zField.modify_on_host();
  }
}

} // mesh
} // stk

#endif // STK_MESH_BASE_NGPFIELDBLAS_HPP

