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
#include "Kokkos_DualView.hpp"
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
    const EXEC_SPACE& /*execSpace*/,
    bool isDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  return operate_on_ngp_mesh<EXEC_SPACE>() || isDeviceExecSpaceUserOverride;
}

template <typename EXEC_SPACE>
void mark_field_modified(const mesh::FieldBase& field, EXEC_SPACE execSpace, bool isDeviceExecSpaceUserOverride)
{
  field.clear_sync_state();
  if (ngp_field_blas::impl::mark_modified_on_device(execSpace, isDeviceExecSpaceUserOverride)) {
    field.modify_on_device();
  }
  else {
    field.modify_on_host();
  }
}

template<typename FieldViewType>
void construct_device_fields(FieldViewType& ngpFields)
{
  using NgpFieldType = typename FieldViewType::value_type;
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, ngpFields.size()),
                       KOKKOS_LAMBDA(const unsigned& i) {
                         new (&ngpFields(i)) NgpFieldType();
                       });
}

template<class Scalar, class NGP_FIELD_TYPE, typename ExecSpace>
class FieldFill {
public:
  FieldFill(const NGP_FIELD_TYPE* fields, int fieldCount, Scalar inputAlpha)
  : ngpFieldsDynamic("ngp_fields", 0), alpha(inputAlpha), nfields(fieldCount)
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        ngpFieldsStatic[i] = fields[i];
      }
    } else
    {
      constexpr bool accessible = Kokkos::SpaceAccessibility<stk::ngp::HostExecSpace, MemorySpace>::accessible;
  
      if (accessible) {
        Kokkos::resize(ngpFieldsDynamic, nfields);
      }
      else {
        Kokkos::resize(Kokkos::WithoutInitializing, ngpFieldsDynamic, nfields);
        construct_device_fields(ngpFieldsDynamic);
      }

      auto ngpFieldsDynamicHost = Kokkos::create_mirror_view(ngpFieldsDynamic);

      for (int i=0; i < nfields; ++i)
      {
        ngpFieldsDynamicHost[i] = fields[i];
      }
      Kokkos::deep_copy(ngpFieldsDynamic, ngpFieldsDynamicHost);
    }
  }

  KOKKOS_FUNCTION ~FieldFill() { }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entityIndex) const
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponents(ngpFieldsStatic[i], entityIndex);
      }
    } else
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponents(ngpFieldsDynamic[i], entityIndex);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void setComponents(const NGP_FIELD_TYPE& ngpField, const stk::mesh::FastMeshIndex& entityIndex) const
  {
    const int numComponents = ngpField.get_num_components_per_entity(entityIndex);
    for(int component=0; component<numComponents; ++component) {
      ngpField(entityIndex, component) = alpha;
    }
  }

  using FieldView = Kokkos::View<NGP_FIELD_TYPE*, ExecSpace>;
  using FieldHostView = typename FieldView::host_mirror_type;
  using MemorySpace = typename FieldView::traits::memory_space;
  static constexpr int STATIC_FIELD_LIMIT = 4;
  NGP_FIELD_TYPE ngpFieldsStatic[STATIC_FIELD_LIMIT];
  FieldView ngpFieldsDynamic;
  Scalar alpha;
  int nfields;
};

template<class Scalar, typename NGP_FIELD_TYPE, typename ExecSpace>
class FieldFillComponent {
public:
  static_assert(std::is_same_v<NGP_FIELD_TYPE, stk::mesh::HostField<Scalar>> ||
                std::is_same_v<NGP_FIELD_TYPE, stk::mesh::DeviceField<Scalar>>);
  FieldFillComponent(const NGP_FIELD_TYPE* fields, int fieldCount, Scalar inputAlpha, int inputComponent)
  : ngpFieldsDynamic("ngp_fields", 0), alpha(inputAlpha), component(inputComponent), nfields(fieldCount)
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        ngpFieldsStatic[i] = fields[i];
      }
    } else
    {
      constexpr bool accessible = Kokkos::SpaceAccessibility<stk::ngp::HostExecSpace, MemorySpace>::accessible;
  
      if constexpr (accessible) {
        Kokkos::resize(ngpFieldsDynamic, nfields);
      }
      else {
        Kokkos::resize(Kokkos::WithoutInitializing, ngpFieldsDynamic, nfields);
        construct_device_fields(ngpFieldsDynamic);
      }

      auto ngpFieldsDynamicHost = Kokkos::create_mirror_view(ngpFieldsDynamic);

      for (int i=0; i < nfields; ++i)
      {
        ngpFieldsDynamicHost(i) = fields[i];
      }

      Kokkos::deep_copy(ngpFieldsDynamic, ngpFieldsDynamicHost);
    }
  }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entityIndex) const
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponent(ngpFieldsStatic[i], entityIndex);
      }
    } else
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponent(ngpFieldsDynamic(i), entityIndex);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void setComponent(const NGP_FIELD_TYPE& ngpField, const stk::mesh::FastMeshIndex& entityIndex) const
  {
      for (int i=0; i < nfields; ++i)
      {
        const int numComponents = ngpField.get_num_components_per_entity(entityIndex);
        STK_NGP_ThrowRequire(component < numComponents);
        ngpField(entityIndex, component) = alpha;
      }
  }

  using FieldView = Kokkos::View<NGP_FIELD_TYPE*, ExecSpace>;
  using FieldHostView = typename FieldView::host_mirror_type;
  using MemorySpace = typename FieldView::traits::memory_space;
  static constexpr int STATIC_FIELD_LIMIT = 4;
  NGP_FIELD_TYPE ngpFieldsStatic[STATIC_FIELD_LIMIT];
  FieldView ngpFieldsDynamic;
  Scalar alpha;
  int component;
  int nfields;
};

template<class Scalar, class NGP_MESH_TYPE, class NGP_FIELD_TYPE, class EXEC_SPACE>
void field_fill_for_each_entity(const NGP_MESH_TYPE& ngpMesh,
                                const NGP_FIELD_TYPE* ngpFields,
                                int nfields,
                                Scalar alpha,
                                int component,
                                const stk::mesh::Selector& selector,
                                const EXEC_SPACE& execSpace)
{
  if (component == -1) {
    FieldFill<Scalar, NGP_FIELD_TYPE, EXEC_SPACE> fieldFill(ngpFields, nfields, alpha);
    stk::mesh::for_each_entity_run(ngpMesh, ngpFields[0].get_rank(), selector, fieldFill, execSpace);
    Kokkos::fence();
  }
  else {
    FieldFillComponent<Scalar, NGP_FIELD_TYPE, EXEC_SPACE> fieldFill(ngpFields, nfields, alpha, component);
    stk::mesh::for_each_entity_run(ngpMesh, ngpFields[0].get_rank(), selector, fieldFill, execSpace);
    Kokkos::fence();
  }
}

template<class Scalar, typename EXEC_SPACE>
void field_fill_impl(const Scalar alpha,
                     const stk::mesh::FieldBase*const* fields,
                     int nfields,
                     int component,
                     const stk::mesh::Selector* selectorPtr,
                     const EXEC_SPACE& execSpace,
                     bool isDeviceExecSpaceUserOverride)
{
  STK_ThrowRequireMsg(nfields > 0, "must have one or more fields");
  for (int i=0; i < nfields; ++i)
  {
    fields[i]->clear_sync_state();
  }

  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(*fields[0]));
    for (int i=1; i < nfields; ++i) {
      *(fieldSelector.get()) &= stk::mesh::Selector(*fields[i]);
    }
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(fields[0]->get_mesh());
    if (nfields == 1)
    {
      stk::mesh::NgpField<Scalar> ngpField = stk::mesh::get_updated_ngp_field<Scalar>(*fields[0]);
      field_fill_for_each_entity(ngpMesh, &ngpField, nfields, alpha, component, selector, execSpace);
    } else
    {
      std::vector<stk::mesh::NgpField<Scalar>> ngpFields;
      for (int i=0; i < nfields; ++i)
      {
        ngpFields.push_back(stk::mesh::get_updated_ngp_field<Scalar>(*fields[i]));
      }
      field_fill_for_each_entity(ngpMesh, ngpFields.data(), nfields, alpha, component, selector, execSpace);
    }
  }
  else {
    if (nfields == 1) {
      stk::mesh::field_fill(alpha, *fields[0], selector);
    }
    else {
      std::vector<const stk::mesh::FieldBase*> fieldsVec(fields, fields+nfields);
      stk::mesh::field_fill(alpha, fieldsVec, selector);
    }
  }

  for (int i=0; i < nfields; ++i)
  {
    mark_field_modified(*fields[i], execSpace, isDeviceExecSpaceUserOverride);
  }
}

template <class Scalar, class NGP_FIELD_TYPE>
class FieldAMax
{
 public:
  FieldAMax(const NGP_FIELD_TYPE& ngpXfield) : ngpX(ngpXfield), resultDevice("deviceResult", 1) {}

  KOKKOS_INLINE_FUNCTION Scalar abs(const Scalar& x) const { return x > 0 ? x : -x; }

  Scalar get_amax_val() const
  {
    auto resultHost = Kokkos::create_mirror_view(resultDevice);
    Kokkos::deep_copy(resultHost, resultDevice);
    return resultHost(0);
  }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entityIndex) const
  {
    const unsigned numComponents = ngpX.get_num_components_per_entity(entityIndex);
    for (unsigned d = 0; d < numComponents; ++d) {
      Scalar fieldVal = ngpX(entityIndex, d);
      Scalar absFieldVal = abs(fieldVal);
      Kokkos::atomic_max(&resultDevice(0), absFieldVal);
    }
  }

  NGP_FIELD_TYPE ngpX;
  Kokkos::View<Scalar*, stk::ngp::ExecSpace> resultDevice;
};

template <class Scalar, class EXEC_SPACE>
Scalar field_amax_no_mark_t(
    const stk::mesh::FieldBase& xField, const stk::mesh::Selector& selector, const EXEC_SPACE& /*execSpace*/)
{
  Scalar amaxOut = 0;

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    xField.sync_to_device();
    stk::mesh::NgpField<Scalar>& ngpX = stk::mesh::get_updated_ngp_field<Scalar>(xField);
    auto ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
    FieldAMax<Scalar, stk::mesh::NgpField<Scalar>> fieldAMax(ngpX);
    stk::mesh::for_each_entity_run(ngpMesh, xField.entity_rank(), selector, fieldAMax);

    Scalar localAmax = fieldAMax.get_amax_val();
    auto comm = xField.get_mesh().parallel();
    stk::all_reduce_max(comm, &localAmax, &amaxOut, 1u);
  } else {
    xField.sync_to_host();
    double tmpAmax = 0.0;
    stk::mesh::field_amax(tmpAmax, xField, selector);
    amaxOut = tmpAmax;
  }

  return amaxOut;
}

template <class Scalar, class EXEC_SPACE>
void field_amax_no_mark_mod(Scalar& amaxOut,
    const stk::mesh::FieldBase& xField,
    const stk::mesh::Selector& selector,
    const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    amaxOut = field_amax_no_mark_t<double>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(float)) {
    amaxOut = field_amax_no_mark_t<float>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(int)) {
    amaxOut = field_amax_no_mark_t<int>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(unsigned)) {
    amaxOut = field_amax_no_mark_t<unsigned>(xField, selector, execSpace);
  } else {
    STK_ThrowErrorMsg("field_amax doesn't yet support fields of type " << dataTraits.type_info.name());
  }
}

template <class Scalar, class EXEC_SPACE>
void field_amax_impl(Scalar& amaxOut,
    const stk::mesh::FieldBase& xField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace,
    bool /*isDeviceExecSpaceUserOverride*/)
{
  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  field_amax_no_mark_mod(amaxOut, xField, selector, execSpace);
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
                          const EXEC_SPACE& /*execSpace*/)
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
    stk::mesh::field_copy(xField, yField, selector);
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

  mark_field_modified(yField, execSpace, isDeviceExecSpaceUserOverride);
}

template <typename Functor>
void apply_functor_on_field(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & zField,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    typename Functor::value_type alpha,
    typename Functor::value_type beta,
    const stk::mesh::Selector & select)
{
  using DataType = typename Functor::value_type;
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

template <typename Scalar>
struct FieldAXPBYZFunctor
{
  using value_type = Scalar;

  FieldAXPBYZFunctor(stk::mesh::NgpMesh & my_ngp_mesh,
      stk::mesh::NgpField<Scalar> & output_field,
      stk::mesh::NgpField<Scalar> & input_x,
      stk::mesh::NgpField<Scalar> & input_y,
      const Scalar & alpha,
      const Scalar & beta)
      : my_mesh(my_ngp_mesh), output(output_field), x(input_x), y(input_y), a(alpha), b(beta)
  {
  }
  stk::mesh::NgpMesh my_mesh;
  stk::mesh::NgpField<Scalar> output;
  stk::mesh::NgpField<Scalar> x;
  stk::mesh::NgpField<Scalar> y;
  const Scalar a;
  const Scalar b;
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

template <typename Scalar>
struct FieldProductFunctor
{
  using value_type = Scalar;

  FieldProductFunctor(stk::mesh::NgpMesh & my_ngp_mesh,
      stk::mesh::NgpField<Scalar> & output_field,
      stk::mesh::NgpField<Scalar> & input_x,
      stk::mesh::NgpField<Scalar> & input_y,
      const Scalar /*a*/,
      const Scalar /*b*/)
      : my_mesh(my_ngp_mesh),
        output(output_field),
        x(input_x), y(input_y)
  {
  }
  stk::mesh::NgpMesh my_mesh;
  stk::mesh::NgpField<Scalar> output;
  stk::mesh::NgpField<Scalar> x;
  stk::mesh::NgpField<Scalar> y;

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
      output.get(f, i) = x.get(f, i) * y.get(f, i);
    }
  }
};

template <typename Scalar>
struct FieldScaleFunctor
{
  using value_type = Scalar;

  FieldScaleFunctor(stk::mesh::NgpMesh & my_ngp_mesh,
      stk::mesh::NgpField<Scalar> & output_field,
      stk::mesh::NgpField<Scalar> & unused_field1,
      stk::mesh::NgpField<Scalar> & unused_field2,
      const Scalar a,
      const Scalar b)
      : my_mesh(my_ngp_mesh),
        output(output_field),
        alpha(a),
        m_unused_field1(unused_field1), m_unused_field2(unused_field2),
        m_beta(b)

  {
  }
  stk::mesh::NgpMesh my_mesh;
  stk::mesh::NgpField<Scalar> output;
  Scalar alpha;

  KOKKOS_FUNCTION
  void operator()(stk::mesh::FastMeshIndex f) const
  {
    unsigned num_components = output.get_num_components_per_entity(f);
    for (unsigned i = 0; i < num_components; ++i)
    {
      output.get(f, i) = alpha * output.get(f, i);
    }
  }

private:
  stk::mesh::NgpField<Scalar> m_unused_field1;
  stk::mesh::NgpField<Scalar> m_unused_field2;
  Scalar m_beta;
};

template <typename Scalar>
struct FieldSwapFunctor
{
  using value_type = Scalar;

  FieldSwapFunctor(stk::mesh::NgpMesh & my_ngp_mesh,
      stk::mesh::NgpField<Scalar> & xFieldInput,
      stk::mesh::NgpField<Scalar> & yFieldInput,
      stk::mesh::NgpField<Scalar> & unused_field2,
      const Scalar a,
      const Scalar b)
      : my_mesh(my_ngp_mesh),
        xField(xFieldInput),
        yField(yFieldInput),
        m_unused_field2(unused_field2),
        m_alpha(a),
        m_beta(b)

  {
  }

  stk::mesh::NgpMesh my_mesh;
  stk::mesh::NgpField<Scalar> xField;
  stk::mesh::NgpField<Scalar> yField;

  KOKKOS_FUNCTION
  void operator()(stk::mesh::FastMeshIndex f) const
  {
    unsigned num_components = xField.get_num_components_per_entity(f);
    unsigned other = yField.get_num_components_per_entity(f);
    num_components = (other < num_components) ? other : num_components;
    for (unsigned i = 0; i < num_components; ++i)
    {
      Scalar tmp = xField.get(f, i);
      xField.get(f, i) = yField.get(f, i);
      yField.get(f, i) = tmp;
    }
  }

private:
  stk::mesh::NgpField<Scalar> m_unused_field2;
  Scalar m_alpha;
  Scalar m_beta;
};

template <typename FieldDataType, typename DataType>
void compute_field_dot(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    DataType& returnVal,
    const stk::mesh::Selector & select)
{
  const stk::mesh::Selector selector = select & stk::mesh::selectField(xField) & stk::mesh::selectField(yField);
  stk::mesh::EntityRank entityRank = xField.entity_rank();
  xField.sync_to_device();
  yField.sync_to_device();
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto& ngpXField = stk::mesh::get_updated_ngp_field<FieldDataType>(xField);
  auto& ngpYField = stk::mesh::get_updated_ngp_field<FieldDataType>(yField);

  using MemSpace = stk::ngp::MemSpace;
  using Layout = Kokkos::LayoutLeft;
  Kokkos::DualView<DataType*, Layout, MemSpace> typeView;
  Kokkos::resize(typeView, 1);
  typeView.modify_device();
  stk::mesh::for_each_entity_run(ngpMesh, entityRank, selector,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
    unsigned num_components = ngpXField.get_num_components_per_entity(entity);
    unsigned other = ngpYField.get_num_components_per_entity(entity);
    num_components = (other < num_components) ? other : num_components;
    for (unsigned i = 0; i < num_components; ++i)
    {
      DataType multVal = ngpXField.get(entity, i) * ngpYField.get(entity, i);
      Kokkos::atomic_add(&typeView.d_view(0), multVal);
    }
  });

  typeView.sync_host();
  returnVal = typeView.h_view(0);

}

template<class DataType, typename EXEC_SPACE>
inline void field_axpbyz_impl(const stk::mesh::BulkData& mesh,
    const DataType alpha,
    const stk::mesh::FieldBase & xField,
    const DataType beta,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::FieldBase & zField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();
  STK_ThrowRequireMsg(dataTraits == yField.data_traits(), "xField and yField must have same datatype");
  STK_ThrowRequireMsg(dataTraits == zField.data_traits(), "xField and zField must have same datatype");

  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField) & stk::mesh::Selector(yField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    if (dataTraits.type_info == typeid(double)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldAXPBYZFunctor<double>>(
        mesh, zField, xField, yField, alpha, beta, selector);
    } else if (dataTraits.type_info == typeid(float)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldAXPBYZFunctor<float>>(
        mesh, zField, xField, yField, alpha, beta, selector);
    } else if (dataTraits.type_info == typeid(int)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldAXPBYZFunctor<int>>(
        mesh, zField, xField, yField, alpha, beta, selector);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldAXPBYZFunctor<unsigned>>(
        mesh, zField, xField, yField, alpha, beta, selector);
    } else
    {
      STK_ThrowErrorMsg("axpbyz doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    xField.sync_to_host();
    yField.sync_to_host();
    stk::mesh::field_copy(yField, zField, selector);
    stk::mesh::field_axpby(alpha, xField, beta, zField, selector);
  }

  mark_field_modified(zField, execSpace, isDeviceExecSpaceUserOverride);
}

template<typename EXEC_SPACE>
inline void field_product_impl(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::FieldBase & zField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();
  STK_ThrowRequireMsg(dataTraits == yField.data_traits(), "xField and yField must have same datatype");
  STK_ThrowRequireMsg(dataTraits == zField.data_traits(), "xField and zField must have same datatype");

  stk::mesh::Selector fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = stk::mesh::Selector(xField) & stk::mesh::Selector(yField);
  } else
  {
    fieldSelector = *selectorPtr;
  }

  if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    if (dataTraits.type_info == typeid(double)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldProductFunctor<double>>(
        mesh, zField, xField, yField, 1, 1, fieldSelector);
    } else if (dataTraits.type_info == typeid(float)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldProductFunctor<float>>(
        mesh, zField, xField, yField, 1, 1, fieldSelector);
    } else if (dataTraits.type_info == typeid(int)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldProductFunctor<int>>(
        mesh, zField, xField, yField, 1, 1, fieldSelector);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldProductFunctor<unsigned>>(
        mesh, zField, xField, yField, 1, 1, fieldSelector);
    } else
    {
      STK_ThrowErrorMsg("field_product doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    xField.sync_to_host();
    yField.sync_to_host();
    stk::mesh::field_product(xField, yField, zField, fieldSelector);
  }

  mark_field_modified(zField, execSpace, isDeviceExecSpaceUserOverride);
}

template<typename Scalar, typename EXEC_SPACE>
inline void field_scale_impl(const stk::mesh::BulkData& mesh,
    const Scalar alpha,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  stk::mesh::Selector fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = stk::mesh::Selector(xField);
  } else
  {
    fieldSelector = *selectorPtr;
  }

  if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    if (dataTraits.type_info == typeid(double)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldScaleFunctor<double>>(
        mesh, xField, xField, xField, alpha, alpha, fieldSelector);
    } else if (dataTraits.type_info == typeid(float)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldScaleFunctor<float>>(
        mesh, xField, xField, xField, alpha, alpha, fieldSelector);
    } else if (dataTraits.type_info == typeid(int)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldScaleFunctor<int>>(
        mesh, xField, xField, xField, alpha, alpha, fieldSelector);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldScaleFunctor<unsigned>>(
        mesh, xField, xField, xField, alpha, alpha, fieldSelector);
    } else
    {
      STK_ThrowErrorMsg("field_product doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    xField.sync_to_host();
    stk::mesh::field_scale(alpha, xField, fieldSelector);
  }

  mark_field_modified(xField, execSpace, isDeviceExecSpaceUserOverride);
}

template<typename EXEC_SPACE>
inline void field_swap_impl(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace,
    bool isDeviceExecSpaceUserOverride)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();
  STK_ThrowRequireMsg(dataTraits == yField.data_traits(), "xField and yField must have same datatype");

  stk::mesh::Selector fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = stk::mesh::Selector(xField);
  } else
  {
    fieldSelector = *selectorPtr;
  }

  if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    if (dataTraits.type_info == typeid(double)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldSwapFunctor<double>>(
        mesh, xField, yField, xField, 0, 0, fieldSelector);
    } else if (dataTraits.type_info == typeid(float)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldSwapFunctor<float>>(
        mesh, xField, yField, xField, 0, 0, fieldSelector);
    } else if (dataTraits.type_info == typeid(int)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldSwapFunctor<int>>(
        mesh, xField, yField, xField, 0, 0, fieldSelector);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      ngp_field_blas::impl::apply_functor_on_field<FieldSwapFunctor<unsigned>>(
        mesh, xField, yField, xField, 0, 0, fieldSelector);
    } else
    {
      STK_ThrowErrorMsg("field_product doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    xField.sync_to_host();
    yField.sync_to_host();
    stk::mesh::field_swap(xField, yField, fieldSelector);
  }

  mark_field_modified(xField, execSpace, isDeviceExecSpaceUserOverride);
  mark_field_modified(yField, execSpace, isDeviceExecSpaceUserOverride);

}

template<class DataType, typename EXEC_SPACE>
inline void field_dot_impl(const stk::mesh::BulkData& mesh,
    const stk::mesh::FieldBase & xField,
    const stk::mesh::FieldBase & yField,
    DataType& returnVal,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();
  STK_ThrowRequireMsg(dataTraits == yField.data_traits(), "xField and yField must have same datatype");

  stk::mesh::Selector fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = stk::mesh::Selector(xField);
  } else
  {
    fieldSelector = *selectorPtr;
  }

  if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    if (dataTraits.type_info == typeid(double)) {
      ngp_field_blas::impl::compute_field_dot<double>(mesh, xField, yField, returnVal, fieldSelector);
    } else if (dataTraits.type_info == typeid(float)) {
      ngp_field_blas::impl::compute_field_dot<float>(mesh, xField, yField, returnVal, fieldSelector);
    } else if (dataTraits.type_info == typeid(int)) {
      ngp_field_blas::impl::compute_field_dot<int>(mesh, xField, yField, returnVal, fieldSelector);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      ngp_field_blas::impl::compute_field_dot<unsigned>(mesh, xField, yField, returnVal, fieldSelector);
    } else
    {
      STK_ThrowErrorMsg("field_dot doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
    DataType globalReturnVal = returnVal;
    stk::all_reduce_sum(mesh.parallel(), &returnVal, &globalReturnVal, 1u);
    returnVal = globalReturnVal;
  }
  else {
    xField.sync_to_host();
    yField.sync_to_host();
    stk::mesh::field_dot(returnVal, xField, yField, *selectorPtr, mesh.parallel());
  }
}


//************ end of implementation detail *********************************

} // namespace impl
} // ngp_field_blas
} // stk

#endif // STK_MESH_BASEIMPL_NGPFIELDBLASIMPL_HPP
