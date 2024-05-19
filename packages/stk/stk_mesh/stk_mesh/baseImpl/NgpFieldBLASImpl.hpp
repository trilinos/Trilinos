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

#ifndef STK_MESH_BASEIMPL_NGPFIELDBLASIMPL_HPP
#define STK_MESH_BASEIMPL_NGPFIELDBLASIMPL_HPP

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/DataTraits.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <complex>
#include <string>
#include <iostream>
#include <algorithm>

namespace stk {
namespace mesh {

namespace impl {
//************ implementation detail, not for public use ********************
//************ public functions are in stk_mesh/base/NgpFieldBLAS.hpp ************

template<class Scalar, class EXEC_SPACE>
inline
void field_fill_no_sync_or_mark(const Scalar alpha, DeviceField<Scalar>& ngpField, const EXEC_SPACE& execSpace)
{
  using FieldExecSpace = typename DeviceField<Scalar>::ExecSpace::execution_space;
  static_assert(Kokkos::SpaceAccessibility<EXEC_SPACE, typename FieldExecSpace::memory_space>::accessible);
  auto ngpView = impl::get_device_data(ngpField);
  Kokkos::deep_copy(execSpace, ngpView, alpha);
}

template<class Scalar, class EXEC_SPACE>
inline
void field_fill_no_sync_or_mark(const Scalar alpha, DeviceField<Scalar>& ngpField, const Selector& selector, const EXEC_SPACE& execSpace)
{
  using FieldExecSpace = typename DeviceField<Scalar>::ExecSpace::execution_space;
  static_assert(Kokkos::SpaceAccessibility<EXEC_SPACE, typename FieldExecSpace::memory_space>::accessible);
  auto ngpMesh = get_updated_ngp_mesh(ngpField.get_field_base()->get_mesh());
  for_each_entity_run(ngpMesh, ngpField.get_rank(), selector, KOKKOS_LAMBDA(const FastMeshIndex& entityIndex) {
    const unsigned numComponents = ngpField.get_num_components_per_entity(entityIndex);
    for(unsigned d=0; d<numComponents; ++d) {
      ngpField(entityIndex, d) = alpha;
    }
  });
}

template<class Scalar, class EXEC_SPACE>
inline
void field_fill_no_sync_or_mark(const Scalar alpha, HostField<Scalar>& ngpField, const EXEC_SPACE& execSpace)
{
  using FieldExecSpace = typename DeviceField<Scalar>::ExecSpace::execution_space;
  static_assert(Kokkos::SpaceAccessibility<EXEC_SPACE, typename FieldExecSpace::memory_space>::accessible);
  stk::mesh::field_fill(alpha, *ngpField.get_field_base());
}

template<class Scalar, class EXEC_SPACE>
inline
void field_fill_no_sync_or_mark(const Scalar alpha, HostField<Scalar>& ngpField, const Selector& selector, const EXEC_SPACE& execSpace)
{
  using FieldExecSpace = typename DeviceField<Scalar>::ExecSpace::execution_space;
  static_assert(Kokkos::SpaceAccessibility<EXEC_SPACE, typename FieldExecSpace::memory_space>::accessible);
  stk::mesh::field_fill(alpha, *ngpField.get_field_base(), selector);
}

template<class Scalar, class EXEC_SPACE>
void field_copy_no_sync_or_mark_t(const FieldBase& xField, const FieldBase& yField, const EXEC_SPACE& execSpace)
{
#ifdef STK_USE_DEVICE_MESH
  NgpField<Scalar>& ngpX = get_updated_ngp_field<Scalar>(xField);
  NgpField<Scalar>& ngpY = get_updated_ngp_field<Scalar>(yField);
  using FieldExecSpace = typename DeviceField<Scalar>::ExecSpace::execution_space;
  static_assert(Kokkos::SpaceAccessibility<EXEC_SPACE, typename FieldExecSpace::memory_space>::accessible);
  auto ngpViewX = impl::get_device_data(ngpX);
  auto ngpViewY = impl::get_device_data(ngpY);
  Kokkos::deep_copy(execSpace, ngpViewY, ngpViewX);
#else
  STK_ThrowErrorMsg("field_copy_no_sync_or_mark_t: there should be no way to get here if STK_USE_DEVICE_MESH not defined");
#endif
}

template<class Scalar, class EXEC_SPACE>
void field_copy_no_sync_or_mark_t(const FieldBase& xField, const FieldBase& yField, const Selector& selector, const EXEC_SPACE& execSpace)
{
#ifdef STK_USE_DEVICE_MESH
  NgpField<Scalar>& ngpX = get_updated_ngp_field<Scalar>(xField);
  NgpField<Scalar>& ngpY = get_updated_ngp_field<Scalar>(yField);
  using FieldExecSpace = typename DeviceField<Scalar>::ExecSpace::execution_space;
  static_assert(Kokkos::SpaceAccessibility<EXEC_SPACE, typename FieldExecSpace::memory_space>::accessible);
  auto ngpMesh = get_updated_ngp_mesh(xField.get_mesh());
  for_each_entity_run(ngpMesh, xField.entity_rank(), selector, KOKKOS_LAMBDA(const FastMeshIndex& entityIndex) {
    const unsigned numComponents = ngpX.get_num_components_per_entity(entityIndex);
    STK_NGP_ThrowAssert(numComponents == ngpY.get_num_components_per_entity(entityIndex));
    for(unsigned d=0; d<numComponents; ++d) {
      ngpY(entityIndex, d) = ngpX(entityIndex, d);
    }
  });
#else
  STK_ThrowErrorMsg("field_copy_no_sync_or_mark_t: there should be no way to get here if STK_USE_DEVICE_MESH not defined");
#endif
}

template<class EXEC_SPACE>
void field_copy_no_sync_or_mark(const FieldBase& xField, const FieldBase& yField, const EXEC_SPACE& execSpace)
{
  const DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    field_copy_no_sync_or_mark_t<double>(xField, yField, execSpace);
  }
  else if (dataTraits.type_info == typeid(float)) {
    field_copy_no_sync_or_mark_t<float>(xField, yField, execSpace);
  }
  else if (dataTraits.type_info == typeid(int)) {
    field_copy_no_sync_or_mark_t<int>(xField, yField, execSpace);
  }
  else if (dataTraits.type_info == typeid(unsigned)) {
    field_copy_no_sync_or_mark_t<unsigned>(xField, yField, execSpace);
  }
  else {
    STK_ThrowErrorMsg("field_copy doesn't yet support fields of type "<<dataTraits.type_info.name());
  }
}

template<class EXEC_SPACE>
void field_copy_no_sync_or_mark(const FieldBase& xField, const FieldBase& yField, const Selector& selector, const EXEC_SPACE& execSpace)
{
  const DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    field_copy_no_sync_or_mark_t<double>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(float)) {
    field_copy_no_sync_or_mark_t<float>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(int)) {
    field_copy_no_sync_or_mark_t<int>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(unsigned)) {
    field_copy_no_sync_or_mark_t<unsigned>(xField, yField, selector, execSpace);
  }
  else {
    STK_ThrowErrorMsg("field_copy doesn't yet support fields of type "<<dataTraits.type_info.name());
  }
}

//************ end of implementation detail *********************************

} // namespace impl
} // mesh
} // stk

#endif // STK_MESH_BASEIMPL_NGPFIELDBLASIMPL_HPP

