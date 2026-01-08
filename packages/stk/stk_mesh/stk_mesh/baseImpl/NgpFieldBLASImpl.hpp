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
#include "stk_mesh/base/NgpReductions.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace stk::ngp_field_blas::impl {

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

template<typename Scalar, typename FieldDataType, typename NgpSpace>
class FieldFill {
public:
  FieldFill(const mesh::FieldBase*const* fields, int fieldCount, Scalar inputAlpha)
  : fieldDataDynamic("fieldDataDynamic", 0), alpha(inputAlpha), nfields(fieldCount)
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        fieldDataStatic[i] = fields[i]->template data<Scalar, mesh::OverwriteAll, NgpSpace>();
      }
    } else
    {
      Kokkos::resize(Kokkos::WithoutInitializing, fieldDataDynamic, nfields);

      auto fieldDataDynamicHost = Kokkos::create_mirror_view(fieldDataDynamic);

      for (int i=0; i < nfields; ++i)
      {
        fieldDataDynamicHost(i) = fields[i]->template data<Scalar, mesh::OverwriteAll, NgpSpace>();
      }
      Kokkos::deep_copy(fieldDataDynamic, fieldDataDynamicHost);
    }
  }

  KOKKOS_FUNCTION ~FieldFill() { }

  KOKKOS_FUNCTION
  void operator()(const mesh::FastMeshIndex& entityIndex) const
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponents(fieldDataStatic[i], entityIndex);
      }
    } else
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponents(fieldDataDynamic[i], entityIndex);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void setComponents(const FieldDataType& fieldData, const mesh::FastMeshIndex& entityIndex) const
  {
    auto fieldValues = fieldData.entity_values(entityIndex);
    for (mesh::ScalarIdx scalar : fieldValues.scalars()) {
      fieldValues(scalar) = alpha;
    }
  }

  using FieldDataView = Kokkos::View<FieldDataType*, typename NgpSpace::mem_space>;
  static constexpr int STATIC_FIELD_LIMIT = 4;
  FieldDataType fieldDataStatic[STATIC_FIELD_LIMIT];
  FieldDataView fieldDataDynamic;
  Scalar alpha;
  int nfields;
};

template<class Scalar, typename FieldDataType, typename NgpSpace>
class FieldFillComponent {
public:
  FieldFillComponent(const mesh::FieldBase*const* fields, int fieldCount, Scalar inputAlpha, int inputComponent)
  : fieldDataDynamic("fieldDataDynamic", 0), alpha(inputAlpha), component(inputComponent), nfields(fieldCount)
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        fieldDataStatic[i] = fields[i]->template data<Scalar, mesh::ReadWrite, NgpSpace>();
      }
    } else
    {
      Kokkos::resize(Kokkos::WithoutInitializing, fieldDataDynamic, nfields);

      auto fieldDataDynamicHost = Kokkos::create_mirror_view(fieldDataDynamic);

      for (int i=0; i < nfields; ++i)
      {
        fieldDataDynamicHost(i) = fields[i]->template data<Scalar, mesh::ReadWrite, NgpSpace>();
      }

      Kokkos::deep_copy(fieldDataDynamic, fieldDataDynamicHost);
    }
  }

  KOKKOS_FUNCTION
  void operator()(const mesh::FastMeshIndex& entityIndex) const
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponent(fieldDataStatic[i], entityIndex);
      }
    } else
    {
      for (int i=0; i < nfields; ++i)
      {
        setComponent(fieldDataDynamic(i), entityIndex);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void setComponent(const FieldDataType& fieldData, const mesh::FastMeshIndex& entityIndex) const
  {
    auto fieldValues = fieldData.entity_values(entityIndex);
    fieldValues(mesh::ScalarIdx(component)) = alpha;
  }

  using FieldDataView = Kokkos::View<FieldDataType*, typename NgpSpace::mem_space>;
  static constexpr int STATIC_FIELD_LIMIT = 4;
  FieldDataType fieldDataStatic[STATIC_FIELD_LIMIT];
  FieldDataView fieldDataDynamic;
  Scalar alpha;
  int component;
  int nfields;
};

template <typename NgpSpace, typename Scalar>
void field_fill_for_each_entity(const mesh::FieldBase*const* fields,
                                int nfields,
                                Scalar alpha,
                                int component,
                                const stk::mesh::Selector& selector,
                                const typename NgpSpace::exec_space& execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(fields[0]->get_mesh());
  if (component == -1) {
    using FieldDataType = decltype(fields[0]->template data<Scalar, mesh::OverwriteAll, NgpSpace>());
    auto fieldFill = FieldFill<Scalar, FieldDataType, NgpSpace>(fields, nfields, alpha);
    for_each_entity_run(ngpMesh, fields[0]->entity_rank(), selector, fieldFill, execSpace);
  }
  else {
    using FieldDataType = decltype(fields[0]->template data<Scalar, mesh::ReadWrite, NgpSpace>());
    auto fieldFill = FieldFillComponent<Scalar, FieldDataType, NgpSpace>(fields, nfields, alpha, component);
    for_each_entity_run(ngpMesh, fields[0]->entity_rank(), selector, fieldFill, execSpace);
  }
  Kokkos::fence();
}

template<typename Scalar, typename EXEC_SPACE>
void field_fill_impl(const Scalar alpha,
                     const stk::mesh::FieldBase*const* fields,
                     int nfields,
                     int component,
                     const stk::mesh::Selector* selectorPtr,
                     const EXEC_SPACE& execSpace,
                     bool /* isDeviceExecSpaceUserOverride */)
{
  STK_ThrowRequireMsg(nfields > 0, "must have one or more fields");

  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(*fields[0]));
    for (int i=1; i < nfields; ++i) {
      *(fieldSelector.get()) &= stk::mesh::Selector(*fields[i]);
    }
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    field_fill_for_each_entity<ngp::DeviceSpace>(fields, nfields, alpha, component, selector, execSpace);
  }
  else {
    if (nfields == 1) {
      if (component == -1) {
        stk::mesh::field_fill(alpha, *fields[0], selector);
      }
      else {
        stk::mesh::field_fill(alpha, *fields[0], component, selector);
      }
    }
    else {
      std::vector<const stk::mesh::FieldBase*> fieldsVec(fields, fields+nfields);
      if (component == -1) {
        stk::mesh::field_fill(alpha, fieldsVec, selector);
      }
      else {
        stk::mesh::field_fill(alpha, fieldsVec, component, selector);
      }
    }
  }
}

template <typename Scalar, typename EXEC_SPACE>
Scalar field_amax_on_device_t(const stk::mesh::FieldBase& xField,
                              const stk::mesh::Selector& selector,
                              const EXEC_SPACE& execSpace)
{
  Scalar amaxOut = 0;

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
    auto fieldData = xField.data<Scalar, mesh::ReadOnly, ngp::DeviceSpace>();

    auto result = Kokkos::View<Scalar, ngp::DeviceSpace::mem_space>("result");
    auto maxer = Kokkos::Max<Scalar, EXEC_SPACE>(result);
    mesh::for_each_entity_reduce(ngpMesh,
                                 xField.entity_rank(),
                                 selector,
                                 maxer,
                                 KOKKOS_LAMBDA(mesh::FastMeshIndex entityIndex,
                                               Scalar& localMax) {
                                   auto fieldValues = fieldData.entity_values(entityIndex);
                                   for (mesh::ScalarIdx scalar : fieldValues.scalars()) {
                                     Scalar fieldValue = fieldValues(scalar);
                                     Scalar absValue = fieldValue > 0 ? fieldValue : -fieldValue;
                                     localMax = Kokkos::max(localMax, absValue);
                                   }
                                 });
    Kokkos::fence();
    auto hostResult = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
    Scalar localAmax = hostResult();
    auto comm = xField.get_mesh().parallel();
    stk::all_reduce_max(comm, &localAmax, &amaxOut, 1u);
  } else {
    double tmpAmax = 0.0;
    stk::mesh::field_amax(tmpAmax, xField, selector);
    amaxOut = tmpAmax;
  }

  return amaxOut;
}

template <typename Scalar, typename EXEC_SPACE>
void field_amax_on_device(Scalar& amaxOut,
                          const stk::mesh::FieldBase& xField,
                          const stk::mesh::Selector& selector,
                          const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    amaxOut = field_amax_on_device_t<double>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(float)) {
    amaxOut = field_amax_on_device_t<float>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(int)) {
    amaxOut = field_amax_on_device_t<int>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(unsigned)) {
    amaxOut = field_amax_on_device_t<unsigned>(xField, selector, execSpace);
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

  field_amax_on_device(amaxOut, xField, selector, execSpace);
}

template<class Scalar, class EXEC_SPACE>
void field_copy_on_device_t(const stk::mesh::FieldBase& xField,
                          const stk::mesh::FieldBase& yField,
                          const stk::mesh::Selector& selector,
                          const EXEC_SPACE& execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto fieldDataX = xField.data<Scalar, mesh::ReadOnly, ngp::DeviceSpace>();
  auto fieldDataY = yField.data<Scalar, mesh::OverwriteAll, ngp::DeviceSpace>();

  stk::mesh::for_each_entity_run(ngpMesh,
                                 xField.entity_rank(),
                                 selector,
                                 KOKKOS_LAMBDA(mesh::FastMeshIndex entityIndex) {
                                   auto fieldValuesX = fieldDataX.entity_values(entityIndex);
                                   auto fieldValuesY = fieldDataY.entity_values(entityIndex);
                                   for (mesh::ScalarIdx scalar : fieldValuesX.scalars()) {
                                     fieldValuesY(scalar) = fieldValuesX(scalar);
                                   }
                                 },
                                 execSpace);
  Kokkos::fence();
}

template<class EXEC_SPACE>
void field_copy_on_device(const stk::mesh::FieldBase& xField,
                            const stk::mesh::FieldBase& yField,
                            const stk::mesh::Selector& selector,
                            const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    field_copy_on_device_t<double>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(float)) {
    field_copy_on_device_t<float>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(int)) {
    field_copy_on_device_t<int>(xField, yField, selector, execSpace);
  }
  else if (dataTraits.type_info == typeid(unsigned)) {
    field_copy_on_device_t<unsigned>(xField, yField, selector, execSpace);
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
                     bool /* isDeviceExecSpaceUserOverride */)
{
  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField) & stk::mesh::Selector(yField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    field_copy_on_device(xField, yField, selector, execSpace);
  }
  else {
    stk::mesh::field_copy(xField, yField, selector);
  }
}

template <typename Scalar, typename NgpSpace>
void apply_axpy_on_field(const mesh::BulkData & mesh,
                           const mesh::FieldBase & xField,
                           const mesh::FieldBase & yField,
                           Scalar a,
                           const stk::mesh::Selector & select,
                           typename NgpSpace::exec_space execSpace)
{
  const stk::mesh::Selector selector = select &
                                       stk::mesh::selectField(xField) &
                                       stk::mesh::selectField(yField);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>();
  auto y = yField.data<Scalar, mesh::ReadWrite, NgpSpace>();

  mesh::for_each_entity_run(ngpMesh,
                            yField.entity_rank(),
                            selector,
                            KOKKOS_LAMBDA(mesh::FastMeshIndex f) {
                              auto xValues = x.entity_values(f);
                              auto yValues = y.entity_values(f);
                              for (mesh::ScalarIdx scalar : yValues.scalars()) {
                                yValues(scalar) += a * xValues(scalar);
                              }
                            },
                            execSpace);
  Kokkos::fence();
}

template<class DataType, typename EXEC_SPACE>
inline void field_axpy_impl(const stk::mesh::BulkData& mesh,
                              const DataType alpha,
                              const stk::mesh::FieldBase & xField,
                              const stk::mesh::FieldBase & yField,
                              const stk::mesh::Selector* selectorPtr,
                              const EXEC_SPACE& execSpace,
                              bool /* isDeviceExecSpaceUserOverride */)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();
  STK_ThrowRequireMsg(dataTraits == yField.data_traits(), "xField and yField must have same datatype");

  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField) & stk::mesh::Selector(yField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    if (dataTraits.type_info == typeid(double)) {
      apply_axpy_on_field<double, ngp::DeviceSpace>(mesh, xField, yField, alpha, selector, execSpace);
    } else if (dataTraits.type_info == typeid(float)) {
      apply_axpy_on_field<float, ngp::DeviceSpace>(mesh, xField, yField, alpha, selector, execSpace);
    } else if (dataTraits.type_info == typeid(int)) {
      apply_axpy_on_field<int, ngp::DeviceSpace>(mesh, xField, yField, alpha, selector, execSpace);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      apply_axpy_on_field<unsigned, ngp::DeviceSpace>(mesh, xField, yField, alpha, selector, execSpace);
    } else
    {
      STK_ThrowErrorMsg("axpy doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    stk::mesh::field_axpy(alpha, xField, yField, selector);
  }
}

template <typename Scalar, typename NgpSpace>
void apply_axpby_on_field(const mesh::BulkData & mesh,
                           const mesh::FieldBase & xField,
                           const mesh::FieldBase & yField,
                           Scalar a,
                           Scalar b,
                           const stk::mesh::Selector & select,
                           typename NgpSpace::exec_space execSpace)
{
  const stk::mesh::Selector selector = select &
                                       stk::mesh::selectField(xField) &
                                       stk::mesh::selectField(yField);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>();
  auto y = yField.data<Scalar, mesh::ReadWrite, NgpSpace>();

  mesh::for_each_entity_run(ngpMesh,
                            yField.entity_rank(),
                            selector,
                            KOKKOS_LAMBDA(mesh::FastMeshIndex f) {
                              auto xValues = x.entity_values(f);
                              auto yValues = y.entity_values(f);
                              for (mesh::ScalarIdx scalar : yValues.scalars()) {
                                yValues(scalar) = a * xValues(scalar) + b * yValues(scalar);
                              }
                            },
                            execSpace);
  Kokkos::fence();
}

template<class DataType, typename EXEC_SPACE>
inline void field_axpby_impl(const stk::mesh::BulkData& mesh,
                              const DataType alpha,
                              const stk::mesh::FieldBase & xField,
                              const DataType beta,
                              const stk::mesh::FieldBase & yField,
                              const stk::mesh::Selector* selectorPtr,
                              const EXEC_SPACE& execSpace,
                              bool /* isDeviceExecSpaceUserOverride */)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();
  STK_ThrowRequireMsg(dataTraits == yField.data_traits(), "xField and yField must have same datatype");

  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField) & stk::mesh::Selector(yField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  if constexpr (ngp_field_blas::impl::operate_on_ngp_mesh<EXEC_SPACE>()) {
    if (dataTraits.type_info == typeid(double)) {
      apply_axpby_on_field<double, ngp::DeviceSpace>(mesh, xField, yField, alpha, beta, selector, execSpace);
    } else if (dataTraits.type_info == typeid(float)) {
      apply_axpby_on_field<float, ngp::DeviceSpace>(mesh, xField, yField, alpha, beta, selector, execSpace);
    } else if (dataTraits.type_info == typeid(int)) {
      apply_axpby_on_field<int, ngp::DeviceSpace>(mesh, xField, yField, alpha, beta, selector, execSpace);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      apply_axpby_on_field<unsigned, ngp::DeviceSpace>(mesh, xField, yField, alpha, beta, selector, execSpace);
    } else
    {
      STK_ThrowErrorMsg("axpby doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    stk::mesh::field_axpby(alpha, xField, beta, yField, selector);
  }
}

template <typename Scalar, typename NgpSpace>
void apply_axpbyz_on_field(const mesh::BulkData & mesh,
                           const mesh::FieldBase & zField,
                           const mesh::FieldBase & xField,
                           const mesh::FieldBase & yField,
                           Scalar a,
                           Scalar b,
                           const stk::mesh::Selector & select,
                           typename NgpSpace::exec_space execSpace)
{
  const stk::mesh::Selector selector = select &
                                       stk::mesh::selectField(zField) &
                                       stk::mesh::selectField(xField) &
                                       stk::mesh::selectField(yField);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto output = zField.data<Scalar, mesh::OverwriteAll, NgpSpace>();
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>();
  auto y = yField.data<Scalar, mesh::ReadOnly, NgpSpace>();

  mesh::for_each_entity_run(ngpMesh,
                            zField.entity_rank(),
                            selector,
                            KOKKOS_LAMBDA(mesh::FastMeshIndex f) {
                              auto outputValues = output.entity_values(f);
                              auto xValues = x.entity_values(f);
                              auto yValues = y.entity_values(f);
                              for (mesh::ScalarIdx scalar : outputValues.scalars()) {
                                outputValues(scalar) = a * xValues(scalar) + b * yValues(scalar);
                              }
                            },
                            execSpace);
  Kokkos::fence();
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
                              bool /* isDeviceExecSpaceUserOverride */)
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
      apply_axpbyz_on_field<double, ngp::DeviceSpace>(mesh, zField, xField, yField, alpha, beta, selector, execSpace);
    } else if (dataTraits.type_info == typeid(float)) {
      apply_axpbyz_on_field<float, ngp::DeviceSpace>(mesh, zField, xField, yField, alpha, beta, selector, execSpace);
    } else if (dataTraits.type_info == typeid(int)) {
      apply_axpbyz_on_field<int, ngp::DeviceSpace>(mesh, zField, xField, yField, alpha, beta, selector, execSpace);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      apply_axpbyz_on_field<unsigned, ngp::DeviceSpace>(mesh, zField, xField, yField, alpha, beta, selector, execSpace);
    } else
    {
      STK_ThrowErrorMsg("axpbyz doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    stk::mesh::field_axpbyz(alpha, xField, beta, yField, zField, selector);
  }
}

template <typename Scalar, typename NgpSpace>
void apply_product_on_field(const mesh::BulkData & mesh,
                            const mesh::FieldBase & zField,
                            const mesh::FieldBase & xField,
                            const mesh::FieldBase & yField,
                            const stk::mesh::Selector & select,
                            typename NgpSpace::exec_space execSpace)
{
  const stk::mesh::Selector selector = select &
                                       stk::mesh::selectField(zField) &
                                       stk::mesh::selectField(xField) &
                                       stk::mesh::selectField(yField);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto output = zField.data<Scalar, mesh::OverwriteAll, NgpSpace>();
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>();
  auto y = yField.data<Scalar, mesh::ReadOnly, NgpSpace>();

  mesh::for_each_entity_run(ngpMesh,
                            zField.entity_rank(),
                            selector,
                            KOKKOS_LAMBDA(mesh::FastMeshIndex f) {
                              auto outputValues = output.entity_values(f);
                              auto xValues = x.entity_values(f);
                              auto yValues = y.entity_values(f);
                              for (mesh::ScalarIdx scalar : outputValues.scalars()) {
                                outputValues(scalar) = xValues(scalar) * yValues(scalar);
                              }
                            },
                            execSpace);
  Kokkos::fence();
}

template<typename EXEC_SPACE>
inline void field_product_impl(const stk::mesh::BulkData& mesh,
                               const stk::mesh::FieldBase & xField,
                               const stk::mesh::FieldBase & yField,
                               const stk::mesh::FieldBase & zField,
                               const stk::mesh::Selector* selectorPtr,
                               const EXEC_SPACE& execSpace,
                               bool /* isDeviceExecSpaceUserOverride */)
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
      apply_product_on_field<double, ngp::DeviceSpace>(mesh, zField, xField, yField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(float)) {
      apply_product_on_field<float, ngp::DeviceSpace>(mesh, zField, xField, yField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(int)) {
      apply_product_on_field<int, ngp::DeviceSpace>(mesh, zField, xField, yField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      apply_product_on_field<unsigned, ngp::DeviceSpace>(mesh, zField, xField, yField, fieldSelector, execSpace);
    } else
    {
      STK_ThrowErrorMsg("field_product doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    stk::mesh::field_product(xField, yField, zField, fieldSelector);
  }
}

template <typename Scalar, typename NgpSpace>
void apply_scale_on_field(const mesh::BulkData & mesh,
                          Scalar alpha,
                          const mesh::FieldBase & xField,
                          const stk::mesh::Selector & select,
                          typename NgpSpace::exec_space execSpace)
{
  const stk::mesh::Selector selector = select &
                                       stk::mesh::selectField(xField);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto fieldData = xField.data<Scalar, mesh::ReadWrite, NgpSpace>();

  mesh::for_each_entity_run(ngpMesh,
                            xField.entity_rank(),
                            selector,
                            KOKKOS_LAMBDA(mesh::FastMeshIndex f) {
                              auto fieldValues = fieldData.entity_values(f);
                              for (mesh::ScalarIdx scalar : fieldValues.scalars()) {
                                fieldValues(scalar) *= alpha;
                              }
                            },
                            execSpace);
  Kokkos::fence();
}

template<typename Scalar, typename EXEC_SPACE>
inline void field_scale_impl(const stk::mesh::BulkData& mesh,
                             const Scalar alpha,
                             const stk::mesh::FieldBase & xField,
                             const stk::mesh::Selector* selectorPtr,
                             const EXEC_SPACE& execSpace,
                             bool /* isDeviceExecSpaceUserOverride */)
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
      apply_scale_on_field<double, ngp::DeviceSpace>(mesh, alpha, xField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(float)) {
      apply_scale_on_field<float, ngp::DeviceSpace>(mesh, alpha, xField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(int)) {
      apply_scale_on_field<int, ngp::DeviceSpace>(mesh, alpha, xField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      apply_scale_on_field<unsigned, ngp::DeviceSpace>(mesh, alpha, xField, fieldSelector, execSpace);
    } else
    {
      STK_ThrowErrorMsg("field_product doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    stk::mesh::field_scale(alpha, xField, fieldSelector);
  }
}

template <typename Scalar, typename NgpSpace>
void apply_swap_on_field(const mesh::BulkData & mesh,
                         const mesh::FieldBase & xField,
                         const mesh::FieldBase & yField,
                         const stk::mesh::Selector & select,
                         typename NgpSpace::exec_space execSpace)
{
  const stk::mesh::Selector selector = select &
                                       stk::mesh::selectField(xField) &
                                       stk::mesh::selectField(yField);
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto xFieldData = xField.data<Scalar, mesh::ReadWrite, NgpSpace>();
  auto yFieldData = yField.data<Scalar, mesh::ReadWrite, NgpSpace>();

  mesh::for_each_entity_run(ngpMesh,
                            xField.entity_rank(),
                            selector,
                            KOKKOS_LAMBDA(stk::mesh::FastMeshIndex f) {
                              auto xValues = xFieldData.entity_values(f);
                              auto yValues = yFieldData.entity_values(f);
                              for (mesh::ScalarIdx scalar : xValues.scalars()) {
                                Scalar tmp = xValues(scalar);
                                xValues(scalar) = yValues(scalar);
                                yValues(scalar) = tmp;
                              }
                            },
                            execSpace);
  Kokkos::fence();
}

template<typename EXEC_SPACE>
inline void field_swap_impl(const stk::mesh::BulkData& mesh,
                            const stk::mesh::FieldBase & xField,
                            const stk::mesh::FieldBase & yField,
                            const stk::mesh::Selector* selectorPtr,
                            const EXEC_SPACE& execSpace,
                            bool /* isDeviceExecSpaceUserOverride */)
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
      apply_swap_on_field<double, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(float)) {
      apply_swap_on_field<float, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(int)) {
      apply_swap_on_field<int, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector, execSpace);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      apply_swap_on_field<int, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector, execSpace);
    } else
    {
      STK_ThrowErrorMsg("field_product doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
  }
  else {
    stk::mesh::field_swap(xField, yField, fieldSelector);
  }
}

template <typename Scalar, typename NgpSpace>
Scalar compute_field_dot(const stk::mesh::BulkData & mesh,
                         const stk::mesh::FieldBase & xField,
                         const stk::mesh::FieldBase & yField,
                         const stk::mesh::Selector & select)
{
  const stk::mesh::Selector selector = select &
                                       mesh::selectField(xField) &
                                       mesh::selectField(yField);
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  auto xFieldData = xField.data<Scalar, mesh::ReadOnly, NgpSpace>();
  auto yFieldData = yField.data<Scalar, mesh::ReadOnly, NgpSpace>();

  auto result = Kokkos::View<Scalar, typename NgpSpace::mem_space>("result");
  auto summer = Kokkos::Sum<Scalar, typename NgpSpace::exec_space>(result);
  stk::mesh::for_each_entity_reduce(ngpMesh,
                                    xField.entity_rank(),
                                    selector,
                                    summer,
                                    KOKKOS_LAMBDA(const mesh::FastMeshIndex& entityIndex,
                                                  Scalar& localSum) {
                                      auto xValues = xFieldData.entity_values(entityIndex);
                                      auto yValues = yFieldData.entity_values(entityIndex);
                                      for (mesh::ScalarIdx scalar : xValues.scalars()) {
                                        localSum += xValues(scalar) * yValues(scalar);
                                      }
                                    });
  auto hostResult = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
  return hostResult();

}

template<class Scalar, typename EXEC_SPACE>
inline void field_dot_impl(const stk::mesh::BulkData & mesh,
                           const stk::mesh::FieldBase & xField,
                           const stk::mesh::FieldBase & yField,
                           Scalar& returnVal,
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
      returnVal = compute_field_dot<double, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector);
    } else if (dataTraits.type_info == typeid(float)) {
      returnVal = compute_field_dot<float, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector);
    } else if (dataTraits.type_info == typeid(int)) {
      returnVal = compute_field_dot<int, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector);
    } else if (dataTraits.type_info == typeid(unsigned)) {
      returnVal = compute_field_dot<unsigned, ngp::DeviceSpace>(mesh, xField, yField, fieldSelector);
    } else
    {
      STK_ThrowErrorMsg("field_dot doesn't yet support fields of type "<<dataTraits.type_info.name());
    }
    Scalar globalReturnVal = returnVal;
    stk::all_reduce_sum(mesh.parallel(), &returnVal, &globalReturnVal, 1u);
    returnVal = globalReturnVal;
  }
  else {
    stk::mesh::field_dot(returnVal, xField, yField, *selectorPtr, mesh.parallel());
  }
}

template <typename Scalar, typename EXEC_SPACE>
Scalar field_amin_on_device_t(const stk::mesh::FieldBase& xField,
                              const stk::mesh::Selector& selector,
                              const EXEC_SPACE& execSpace)
{
  Scalar aminOut = 0;

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
    auto fieldData = xField.data<Scalar, mesh::ReadOnly, ngp::DeviceSpace>();

    auto result = Kokkos::View<Scalar, ngp::DeviceSpace::mem_space>("result");
    auto minner = Kokkos::Min<Scalar, EXEC_SPACE>(result);
    mesh::for_each_entity_reduce(ngpMesh,
                                 xField.entity_rank(),
                                 selector,
                                 minner,
                                 KOKKOS_LAMBDA(mesh::FastMeshIndex entityIndex,
                                               Scalar& localMin) {
                                   auto fieldValues = fieldData.entity_values(entityIndex);
                                   for (mesh::ScalarIdx scalar : fieldValues.scalars()) {
                                     Scalar fieldValue = fieldValues(scalar);
                                     Scalar absValue = fieldValue > 0 ? fieldValue : -fieldValue;
                                     localMin = Kokkos::min(localMin, absValue);
                                   }
                                 });
    Kokkos::fence();
    auto hostResult = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
    Scalar localAmin = hostResult();
    auto comm = xField.get_mesh().parallel();
    stk::all_reduce_min(comm, &localAmin, &aminOut, 1u);
  } else {
    double tmpAmin = 0.0;
    stk::mesh::field_amin(tmpAmin, xField, selector);
    aminOut = tmpAmin;
  }

  return aminOut;
}

template <typename Scalar, typename EXEC_SPACE>
void field_amin_on_device(Scalar& aminOut,
                          const stk::mesh::FieldBase& xField,
                          const stk::mesh::Selector& selector,
                          const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    aminOut = field_amin_on_device_t<double>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(float)) {
    aminOut = field_amin_on_device_t<float>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(int)) {
    aminOut = field_amin_on_device_t<int>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(unsigned)) {
    aminOut = field_amin_on_device_t<unsigned>(xField, selector, execSpace);
  } else {
    STK_ThrowErrorMsg("field_amin doesn't yet support fields of type " << dataTraits.type_info.name());
  }
}

template <class Scalar, class EXEC_SPACE>
void field_amin_impl(Scalar& aminOut,
    const stk::mesh::FieldBase& xField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace)
{
  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  field_amin_on_device(aminOut, xField, selector, execSpace);
}

template <typename Scalar, typename EXEC_SPACE>
Scalar field_asum_on_device_t(const stk::mesh::FieldBase& xField,
                               const stk::mesh::Selector& selector,
                               const EXEC_SPACE& execSpace)
{
  Scalar asumOut = 0;

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
    auto fieldData = xField.data<Scalar, mesh::ReadOnly, ngp::DeviceSpace>();

    auto result = Kokkos::View<Scalar, ngp::DeviceSpace::mem_space>("result");
    auto summer = Kokkos::Sum<Scalar, EXEC_SPACE>(result);
    mesh::for_each_entity_reduce(ngpMesh,
                                 xField.entity_rank(),
                                 selector,
                                 summer,
                                 KOKKOS_LAMBDA(mesh::FastMeshIndex entityIndex,
                                               Scalar& localSum) {
                                   auto fieldValues = fieldData.entity_values(entityIndex);
                                   for (mesh::ScalarIdx scalar : fieldValues.scalars()) {
                                     Scalar fieldValue = fieldValues(scalar);
                                     Scalar absValue = fieldValue > 0 ? fieldValue : -fieldValue;
                                     localSum += absValue;
                                   }
                                 });
    Kokkos::fence();
    auto hostResult = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
    Scalar localAsum = hostResult();
    auto comm = xField.get_mesh().parallel();
    stk::all_reduce_sum(comm, &localAsum, &asumOut, 1u);
  } else {
    double tmpAsum = 0.0;
    stk::mesh::field_asum(tmpAsum, xField, selector);
    asumOut = tmpAsum;
  }

  return asumOut;
}

template <typename Scalar, typename EXEC_SPACE>
void field_asum_on_device(Scalar& aminOut,
                          const stk::mesh::FieldBase& xField,
                          const stk::mesh::Selector& selector,
                          const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    aminOut = field_asum_on_device_t<double>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(float)) {
    aminOut = field_asum_on_device_t<float>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(int)) {
    aminOut = field_asum_on_device_t<int>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(unsigned)) {
    aminOut = field_asum_on_device_t<unsigned>(xField, selector, execSpace);
  } else {
    STK_ThrowErrorMsg("field_asum doesn't yet support fields of type " << dataTraits.type_info.name());
  }
}

template <class Scalar, class EXEC_SPACE>
void field_asum_impl(Scalar& asumOut,
    const stk::mesh::FieldBase& xField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace)
{
  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  field_asum_on_device(asumOut, xField, selector, execSpace);
}

template <typename Scalar, typename EXEC_SPACE>
Scalar field_nrm2_on_device_t(const stk::mesh::FieldBase& xField,
                               const stk::mesh::Selector& selector,
                               const EXEC_SPACE& execSpace)
{
  Scalar nrm2Out = 0;

  if constexpr (operate_on_ngp_mesh<EXEC_SPACE>()) {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
    auto fieldData = xField.data<Scalar, mesh::ReadOnly, ngp::DeviceSpace>();

    auto result = Kokkos::View<Scalar, ngp::DeviceSpace::mem_space>("result");
    auto summer = Kokkos::Sum<Scalar, EXEC_SPACE>(result);
    mesh::for_each_entity_reduce(ngpMesh,
                                 xField.entity_rank(),
                                 selector,
                                 summer,
                                 KOKKOS_LAMBDA(mesh::FastMeshIndex entityIndex,
                                               Scalar& localNrm2) {
                                   auto fieldValues = fieldData.entity_values(entityIndex);
                                   for (mesh::ScalarIdx scalar : fieldValues.scalars()) {
                                     Scalar fieldValue = fieldValues(scalar);
                                     localNrm2 += fieldValue * fieldValue;
                                   }
                                 });
    Kokkos::fence();
    auto hostResult = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
    Scalar localNrm2 = hostResult();
    auto comm = xField.get_mesh().parallel();
    stk::all_reduce_sum(comm, &localNrm2, &nrm2Out, 1u);
    nrm2Out = std::sqrt(nrm2Out);
  } else {
    double tmpNrm2 = 0.0;
    stk::mesh::field_nrm2(tmpNrm2, xField, selector);
    nrm2Out = tmpNrm2;
  }

  return nrm2Out;
}

template <typename Scalar, typename EXEC_SPACE>
void field_nrm2_on_device(Scalar& nrm2Out,
                          const stk::mesh::FieldBase& xField,
                          const stk::mesh::Selector& selector,
                          const EXEC_SPACE& execSpace)
{
  const stk::mesh::DataTraits& dataTraits = xField.data_traits();

  if (dataTraits.type_info == typeid(double)) {
    nrm2Out = field_nrm2_on_device_t<double>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(float)) {
    nrm2Out = field_nrm2_on_device_t<float>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(int)) {
    nrm2Out = field_nrm2_on_device_t<int>(xField, selector, execSpace);
  } else if (dataTraits.type_info == typeid(unsigned)) {
    nrm2Out = field_nrm2_on_device_t<unsigned>(xField, selector, execSpace);
  } else {
    STK_ThrowErrorMsg("field_nrm2 doesn't yet support fields of type " << dataTraits.type_info.name());
  }
}

template <class Scalar, class EXEC_SPACE>
void field_nrm2_impl(Scalar& nrm2Out,
    const stk::mesh::FieldBase& xField,
    const stk::mesh::Selector* selectorPtr,
    const EXEC_SPACE& execSpace)
{
  std::unique_ptr<stk::mesh::Selector> fieldSelector;
  if (selectorPtr == nullptr) {
    fieldSelector = std::make_unique<stk::mesh::Selector>(stk::mesh::Selector(xField));
  }
  const stk::mesh::Selector& selector = selectorPtr != nullptr ? *selectorPtr : *(fieldSelector.get());

  field_nrm2_on_device(nrm2Out, xField, selector, execSpace);
}

//************ end of implementation detail *********************************

} // namespace stk::ngp_field_blas::impl

#endif // STK_MESH_BASEIMPL_NGPFIELDBLASIMPL_HPP
