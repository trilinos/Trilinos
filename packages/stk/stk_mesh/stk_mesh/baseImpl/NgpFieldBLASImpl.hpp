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
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <complex>
#include <string>
#include <iostream>
#include <algorithm>

namespace stk {
namespace ngp_field_blas {
namespace impl {

//************ implementation detail, not for public use ********************
//************ public functions are in stk_mesh/base/NgpFieldBLAS.hpp ************

template<class EXEC_SPACE>
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

template<class EXEC_SPACE>
bool mark_modified_on_device(
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  return operate_on_ngp_mesh<EXEC_SPACE>() || isDeviceExecSpaceUserOverride;
}

template<class Scalar, class NGP_FIELD_TYPE>
class FieldFill {
public:
  FieldFill(const NGP_FIELD_TYPE& field, Scalar inputAlpha)
  : ngpField(field), alpha(inputAlpha)
  {}

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entityIndex) const
  {
    const int numComponents = ngpField.get_num_components_per_entity(entityIndex);
    for(int component=0; component<numComponents; ++component) {
      ngpField(entityIndex, component) = alpha;
    }
  }

  NGP_FIELD_TYPE ngpField;
  Scalar alpha;
};

template<class Scalar, class NGP_FIELD_TYPE>
class FieldFillComponent {
public:
  FieldFillComponent(const NGP_FIELD_TYPE& field, Scalar inputAlpha, int inputComponent)
  : ngpField(field), alpha(inputAlpha), component(inputComponent)
  {}

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entityIndex) const
  {
    const int numComponents = ngpField.get_num_components_per_entity(entityIndex);
    STK_NGP_ThrowRequire(component < numComponents);
    ngpField(entityIndex, component) = alpha;
  }

  NGP_FIELD_TYPE ngpField;
  Scalar alpha;
  int component;
};

template<class Scalar, class NGP_MESH_TYPE, class NGP_FIELD_TYPE, class EXEC_SPACE>
void field_fill_for_each_entity(const NGP_MESH_TYPE& ngpMesh,
                                const NGP_FIELD_TYPE& ngpField,
                                Scalar alpha,
                                int component,
                                const stk::mesh::Selector& selector,
                                const EXEC_SPACE& execSpace)
{
  if (component == -1) {
    FieldFill<Scalar, NGP_FIELD_TYPE> fieldFill(ngpField, alpha);
    stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), selector, fieldFill, execSpace);
  }
  else {
    FieldFillComponent<Scalar, NGP_FIELD_TYPE> fieldFill(ngpField, alpha, component);
    stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), selector, fieldFill, execSpace);
  }
}

template<class Scalar, typename EXEC_SPACE>
void field_fill_impl(const Scalar alpha,
                     const stk::mesh::FieldBase& field,
                     int component,
                     const stk::mesh::Selector* selectorPtr,
                     const EXEC_SPACE& execSpace,
                     bool isDeviceExecSpaceUserOverride)
{
  field.clear_sync_state();

  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) { 
    fieldSelector = std::make_unique<stk::mesh::Selector>(field);
  } 
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(field.get_mesh());
    stk::mesh::NgpField<Scalar>& ngpField = stk::mesh::get_updated_ngp_field<Scalar>(field);
    field_fill_for_each_entity(ngpMesh, ngpField, alpha, component, selector, execSpace);
  }
  else {
    stk::mesh::HostMesh hostMesh(field.get_mesh());
    stk::mesh::HostField<Scalar> hostField(field.get_mesh(), field);
    field_fill_for_each_entity(hostMesh, hostField, alpha, component, selector, execSpace);
  }

  if (mark_modified_on_device(execSpace, isDeviceExecSpaceUserOverride)) {
    field.modify_on_device();
  }
  else {
    field.modify_on_host();
  }
}

template<class NGP_FIELD_TYPE>
class FieldCopy {
public:
  FieldCopy(const NGP_FIELD_TYPE& ngpXfield, const NGP_FIELD_TYPE& ngpYfield)
  : ngpX(ngpXfield), ngpY(ngpYfield)
  {}

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entityIndex) const
  {
    const unsigned numComponents = ngpX.get_num_components_per_entity(entityIndex);
    STK_NGP_ThrowAssert(numComponents == ngpY.get_num_components_per_entity(entityIndex));
    for(unsigned d=0; d<numComponents; ++d) {
      ngpY(entityIndex, d) = ngpX(entityIndex, d);
    }
  }

  NGP_FIELD_TYPE ngpX;
  NGP_FIELD_TYPE ngpY;
};

template<class Scalar, class EXEC_SPACE>
void field_copy_no_mark_t(const stk::mesh::FieldBase& xField,
                          const stk::mesh::FieldBase& yField,
                          const stk::mesh::Selector& selector,
                          const EXEC_SPACE& execSpace)
{
  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    xField.sync_to_device();
    stk::mesh::NgpField<Scalar>& ngpX = stk::mesh::get_updated_ngp_field<Scalar>(xField);
    stk::mesh::NgpField<Scalar>& ngpY = stk::mesh::get_updated_ngp_field<Scalar>(yField);
    auto ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
    FieldCopy<stk::mesh::NgpField<Scalar>> fieldCopy(ngpX, ngpY);
    stk::mesh::for_each_entity_run(ngpMesh, xField.entity_rank(), selector, fieldCopy);
  }
  else {
    xField.sync_to_host();
    stk::mesh::HostField<Scalar> hostX(xField.get_mesh(), xField);
    stk::mesh::HostField<Scalar> hostY(yField.get_mesh(), yField);
    stk::mesh::HostMesh hostMesh(xField.get_mesh());
    FieldCopy<stk::mesh::HostField<Scalar>> fieldCopy(hostX, hostY);
    stk::mesh::for_each_entity_run(hostMesh, xField.entity_rank(), selector, fieldCopy);
  }
}

template<class EXEC_SPACE>
void field_copy_no_mark_mod(const stk::mesh::FieldBase& xField,
                            const stk::mesh::FieldBase& yField,
                            const stk::mesh::Selector& selector,
                            const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    field_copy_no_mark_t<double>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(float)) {
    field_copy_no_mark_t<float>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(int)) {
    field_copy_no_mark_t<int>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(unsigned)) {
    field_copy_no_mark_t<unsigned>(xField, yField, selector, execSpace);
  }
  else {
    STK_ThrowErrorMsg("field_copy doesn't yet support fields of type "<<dataTraits.type_info.name());
  }
}

template<class EXEC_SPACE>
void field_copy_impl(const stk::mesh::FieldBase& xField,
                     const stk::mesh::FieldBase& yField,
                     const stk::mesh::Selector* selectorPtr,
                     const EXEC_SPACE& execSpace,
                     bool isDeviceExecSpaceUserOverride)
{
  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) { 
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField) & stk::mesh::Selector(yField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  field_copy_no_mark_mod(xField, yField, selector, execSpace);

  yField.clear_sync_state();
  if (mark_modified_on_device(execSpace, isDeviceExecSpaceUserOverride)) {
    yField.modify_on_device();
  }
  else {
    yField.modify_on_host();
  }
}

template <typename DataType, typename Functor>
void apply_functor_on_field(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & zField,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const DataType alpha,
    const DataType beta,
    const stk::mesh::Selector & select)
{
  const stk::mesh::Selector selector = select & stk::mesh::selectField(zField) &
      stk::mesh::selectField(xField) & stk::mesh::selectField(yField);
  stk::mesh::EntityRank entityRank = zField.entity_rank();
  xField.sync_to_device();
  yField.sync_to_device();
  zField.sync_to_device();
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto& ngpXField = stk::mesh::get_updated_ngp_field<DataType>(xField);
  auto& ngpYField = stk::mesh::get_updated_ngp_field<DataType>(yField);
  auto& ngpFieldResult = stk::mesh::get_updated_ngp_field<DataType>(zField);
  Functor f(ngpMesh, ngpFieldResult, ngpXField, ngpYField, alpha, beta);
  stk::mesh::for_each_entity_run(ngpMesh, entityRank, selector, f);
  zField.modify_on_device();
}


struct FieldAXPBYFunctor
{
  FieldAXPBYFunctor(stk::mesh::NgpMesh & my_ngp_mesh,
      stk::mesh::NgpField<double> & output_field,
      stk::mesh::NgpField<double> & input_x,
      stk::mesh::NgpField<double> & input_y,
      const double & alpha,
      const double & beta)
      : my_mesh(my_ngp_mesh), output(output_field), x(input_x), y(input_y), a(alpha), b(beta)
  {
  }
  stk::mesh::NgpMesh my_mesh;
  stk::mesh::NgpField<double> output;
  stk::mesh::NgpField<double> x;
  stk::mesh::NgpField<double> y;
  const double a;
  const double b;
  KOKKOS_FUNCTION
  void operator()(stk::mesh::FastMeshIndex f) const
  {
    unsigned num_components = output.get_num_components_per_entity(f);
    unsigned other = x.get_num_components_per_entity(f);
    num_components = (other < num_components) ? other : num_components;
    other = y.get_num_components_per_entity(f);
    num_components = (other < num_components) ? other : num_components;
    for (unsigned i = 0; i < num_components; ++i)
    {
      output.get(f, i) = a * x.get(f, i) + b * y.get(f, i);
    }
  }
};

//************ end of implementation detail *********************************

} // namespace impl
} // ngp_field_blas
} // stk

#endif // STK_MESH_BASEIMPL_NGPFIELDBLASIMPL_HPP

