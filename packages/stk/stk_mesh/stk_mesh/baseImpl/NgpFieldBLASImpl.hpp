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
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include "stk_mesh/base/NgpReductions.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include <stk_util/parallel/ParallelReduce.hpp>
#include <Kokkos_StdAlgorithms.hpp>

namespace stk::ngp_field_blas::impl {

//************ implementation detail, not for public use ********************
//************ public functions are in stk_mesh/base/NgpFieldBLAS.hpp ************

template<typename Scalar, typename FieldDataType, typename NgpSpace>
class FieldFill {
public:
  FieldFill(const std::vector<const mesh::FieldBase*>& fields, const stk::NgpVector<unsigned>& inputIds, Scalar inputAlpha, const typename NgpSpace::exec_space& execSpace)
  : bucketIds(inputIds),
    alpha(inputAlpha),
    nfields(fields.size())
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        fieldDataStatic[i] = fields[i]->template data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace);
      }
    } else
    {
      fieldDataDynamic = FieldDataView("fieldDataDynamic", fields.size());
      auto fieldDataDynamicHost = Kokkos::create_mirror_view(stk::ngp::HostPinnedSpace{}, fieldDataDynamic);
      for (int i=0; i < nfields; ++i)
      {
        fieldDataDynamicHost(i) = fields[i]->template data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace);
      }

      Kokkos::Experimental::copy(execSpace, fieldDataDynamicHost, fieldDataDynamic);
    }
  }

  KOKKOS_FUNCTION ~FieldFill() { }

  KOKKOS_FUNCTION
  void operator()(const typename stk::ngp::TeamPolicy<typename NgpSpace::exec_space>::member_type& team) const
  {
    const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank() / nfields);
    const int i = team.league_rank() % nfields;
    auto xValues = (nfields <= STATIC_FIELD_LIMIT) ? fieldDataStatic[i].bucket_values(bucketIndex) : fieldDataDynamic[i].bucket_values(bucketIndex);
    const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
    const unsigned numEntities = xValues.num_entities();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                         [&](const unsigned idx) {
                           const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                           const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                           xValues(entityIdx, scalarIdx) = alpha;
                         });
  }

  using FieldDataView = Kokkos::View<FieldDataType*, typename NgpSpace::mem_space>;
  using HostPinnedFieldDataView = Kokkos::View<FieldDataType*, stk::ngp::HostPinnedSpace>;
  static constexpr int STATIC_FIELD_LIMIT = 4;
  FieldDataType fieldDataStatic[STATIC_FIELD_LIMIT];
  FieldDataView fieldDataDynamic;
  stk::NgpVector<unsigned> bucketIds;
  Scalar alpha;
  int nfields;
};

template<class Scalar, typename FieldDataType, typename NgpSpace>
class FieldFillComponent {
public:
  FieldFillComponent(const std::vector<const mesh::FieldBase*>& fields, const stk::NgpVector<unsigned>& inputIds, Scalar inputAlpha, int inputComponent, const typename NgpSpace::exec_space& execSpace)
  : fieldDataDynamic("fieldDataDynamic", (fields.size() <= STATIC_FIELD_LIMIT) ? 0 : fields.size()),
    bucketIds(inputIds),
    alpha(inputAlpha),
    component(inputComponent),
    nfields(fields.size())
  {
    if (nfields <= STATIC_FIELD_LIMIT)
    {
      for (int i=0; i < nfields; ++i)
      {
        fieldDataStatic[i] = fields[i]->template data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace);
      }
    } else
    {
      auto fieldDataDynamicHost = HostPinnedFieldDataView("fieldDataDynamicHost", fieldDataDynamic.extent(0));
      for (int i=0; i < nfields; ++i)
      {
        fieldDataDynamicHost(i) = fields[i]->template data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace);
      }

      Kokkos::Experimental::copy(execSpace, fieldDataDynamicHost, fieldDataDynamic);
    }
  }

  KOKKOS_FUNCTION
  void operator()(const typename stk::ngp::TeamPolicy<typename NgpSpace::exec_space>::member_type& team) const
  {
    const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank() / nfields);
    const int i = team.league_rank() % nfields;
    auto xValues = (nfields <= STATIC_FIELD_LIMIT) ? fieldDataStatic[i].bucket_values(bucketIndex) : fieldDataDynamic[i].bucket_values(bucketIndex);
    const unsigned numEntities = xValues.num_entities();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numEntities),
                         [&](const unsigned idx) {
                           const auto entityIdx = mesh::EntityIdx(idx);
                           const auto scalarIdx = mesh::ScalarIdx(component);
                           xValues(entityIdx, scalarIdx) = alpha;
                         });
  }

  using FieldDataView = Kokkos::View<FieldDataType*, typename NgpSpace::mem_space>;
  using HostPinnedFieldDataView = Kokkos::View<FieldDataType*, stk::ngp::HostPinnedSpace>;
  static constexpr int STATIC_FIELD_LIMIT = 4;
  FieldDataType fieldDataStatic[STATIC_FIELD_LIMIT];
  FieldDataView fieldDataDynamic;
  stk::NgpVector<unsigned> bucketIds;
  Scalar alpha;
  int component;
  int nfields;
};

template<typename NgpSpace, typename Scalar>
void field_fill_impl(const Scalar alpha,
                     const stk::mesh::FieldBase & xField,
                     int component,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto xData = xField.data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_fill",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         auto xValues = xData.bucket_values(bucketIndex);
                         const unsigned numEntities = xValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numEntities),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx);
                                                const auto scalarIdx = mesh::ScalarIdx(component);
                                                xValues(entityIdx, scalarIdx) = alpha;
                                              });
                       });
}

template<typename NgpSpace, typename Scalar>
void field_fill_impl(const Scalar alpha,
                     const stk::mesh::FieldBase & xField,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto xData = xField.data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_fill",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         auto xValues = xData.bucket_values(bucketIndex);
                         const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
                         const unsigned numEntities = xValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                xValues(entityIdx, scalarIdx) = alpha;
                                              });
                       });
}

template<typename NgpSpace, typename Scalar>
void field_fill_impl(const Scalar alpha,
                     const std::vector<const stk::mesh::FieldBase*> & fields,
                     int component,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(fields.front()->get_mesh());
  using FieldDataType = decltype(fields.front()->template data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace));

  const stk::topology::rank_t rank = fields.front()->entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  auto fieldFill = FieldFillComponent<Scalar, FieldDataType, NgpSpace>(fields, bucketIds, alpha, component, execSpace);
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  Kokkos::parallel_for("field_fill",
                       TeamPolicyType(execSpace, numBuckets * fields.size(), Kokkos::AUTO),
                       fieldFill);
}

template<typename NgpSpace, typename Scalar>
void field_fill_impl(const Scalar alpha,
                     const std::vector<const stk::mesh::FieldBase*> & fields,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(fields.front()->get_mesh());
  using FieldDataType = decltype(fields.front()->template data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace));

  const stk::topology::rank_t rank = fields.front()->entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  auto fieldFill = FieldFill<Scalar, FieldDataType, NgpSpace>(fields, bucketIds, alpha, execSpace);
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  Kokkos::parallel_for("field_fill",
                       TeamPolicyType(execSpace, numBuckets * fields.size(), Kokkos::AUTO),
                       fieldFill);
}

template <typename NgpSpace, typename Scalar>
void field_amax_impl(Scalar& amaxOut,
                     const stk::mesh::FieldBase& xField,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto fieldData = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);

  auto result = Kokkos::View<Scalar, typename NgpSpace::mem_space>("result");
  auto maxer = Kokkos::Max<Scalar, typename NgpSpace::exec_space>(result);
  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  const unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_reduce("field_amax",
                          TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                          KOKKOS_LAMBDA(const TeamHandleType& team, Scalar& teamReduction) {
                            const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                            auto fieldValues = fieldData.bucket_values(bucketIndex);
                            const unsigned numValues = fieldValues.num_entities() * fieldValues.num_scalars();
                            const unsigned numEntities = fieldValues.num_entities();

                            Scalar localReduction;
                            auto threadMaxer = Kokkos::Max<Scalar, typename NgpSpace::exec_space>(localReduction);
                            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numValues),
                                                    [&](const unsigned idx, Scalar& threadReduction) {
                                                      const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                      const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                      Scalar fieldValue = fieldValues(entityIdx, scalarIdx);
                                                      Scalar absValue = fieldValue > 0 ? fieldValue : -fieldValue;
                                                      threadReduction = Kokkos::max(threadReduction, absValue);
                                                    },
                                                    threadMaxer);
                            Kokkos::single(Kokkos::PerTeam(team),
                                           [&]() {
                                             teamReduction = Kokkos::max(teamReduction, localReduction);
                                           });

                          },
                          maxer);
  Kokkos::deep_copy(execSpace, amaxOut, result);
  execSpace.fence();

  Scalar localAmax = amaxOut;
  stk::all_reduce_max(xField.get_mesh().parallel(), &localAmax, &amaxOut, 1u);
}

template<typename NgpSpace, typename Scalar>
void field_copy_on_device(const stk::mesh::FieldBase & xField,
                            const stk::mesh::FieldBase & yField,
                            const stk::mesh::Selector & selector,
                            const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto fieldDataX = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);
  auto fieldDataY = yField.data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_product",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         const auto fieldValuesX = fieldDataX.bucket_values(bucketIndex);
                         auto fieldValuesY = fieldDataY.bucket_values(bucketIndex);
                         const unsigned numValues = fieldValuesX.num_entities() * fieldValuesX.num_scalars();
                         const unsigned numEntities = fieldValuesX.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                fieldValuesY(entityIdx, scalarIdx) = fieldValuesX(entityIdx, scalarIdx);
                                              });
                       });
}

template<typename NgpSpace>
void field_copy_impl(const stk::mesh::FieldBase & xField,
                     const stk::mesh::FieldBase & yField,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  field_datatype_execute(xField,
    [&]<typename T>(const stk::mesh::FieldBase&) {
      field_copy_on_device<NgpSpace, T>(xField, yField, selector, execSpace);
    }
  );
}

template<typename NgpSpace, typename Scalar>
inline void field_axpy_impl(const Scalar alpha,
                            const stk::mesh::FieldBase & xField,
                            const stk::mesh::FieldBase & yField,
                            const stk::mesh::Selector & selector,
                            const typename NgpSpace::exec_space& execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);
  auto y = yField.data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_axpy",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         const auto xValues = x.bucket_values(bucketIndex);
                         auto yValues = y.bucket_values(bucketIndex);
                         const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
                         const unsigned numEntities = xValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                yValues(entityIdx, scalarIdx) += alpha * xValues(entityIdx, scalarIdx);
                                              });
                       });
}

template<typename NgpSpace, typename Scalar>
inline void field_axpby_impl(const Scalar alpha,
                             const stk::mesh::FieldBase & xField,
                             const Scalar beta,
                             const stk::mesh::FieldBase & yField,
                             const stk::mesh::Selector & selector,
                             const typename NgpSpace::exec_space& execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);
  auto y = yField.data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_axpby",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         const auto xValues = x.bucket_values(bucketIndex);
                         auto yValues = y.bucket_values(bucketIndex);
                         const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
                         const unsigned numEntities = xValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                yValues(entityIdx, scalarIdx) = alpha * xValues(entityIdx, scalarIdx) + beta * yValues(entityIdx, scalarIdx);
                                              });
                       });
}

template<typename NgpSpace, typename Scalar>
inline void field_axpbyz_impl(const Scalar alpha,
                              const stk::mesh::FieldBase & xField,
                              const Scalar beta,
                              const stk::mesh::FieldBase & yField,
                              const stk::mesh::FieldBase & zField,
                              const stk::mesh::Selector & selector,
                              const typename NgpSpace::exec_space& execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto output = zField.data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace);
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);
  auto y = yField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_axpbyz",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         auto outputValues = output.bucket_values(bucketIndex);
                         const auto xValues = x.bucket_values(bucketIndex);
                         const auto yValues = y.bucket_values(bucketIndex);
                         const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
                         const unsigned numEntities = xValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                outputValues(entityIdx, scalarIdx) = alpha * xValues(entityIdx, scalarIdx) + beta * yValues(entityIdx, scalarIdx);
                                              });
                       });
}

template <typename NgpSpace, typename Scalar>
void apply_product_on_field(const mesh::FieldBase & xField,
                            const mesh::FieldBase & yField,
                            const mesh::FieldBase & zField,
                            const stk::mesh::Selector & selector,
                            const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto output = zField.data<Scalar, mesh::OverwriteAll, NgpSpace>(execSpace);
  auto x = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);
  auto y = yField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_product",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         auto outputValues = output.bucket_values(bucketIndex);
                         const auto xValues = x.bucket_values(bucketIndex);
                         const auto yValues = y.bucket_values(bucketIndex);
                         const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
                         const unsigned numEntities = xValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                outputValues(entityIdx, scalarIdx) = xValues(entityIdx, scalarIdx) * yValues(entityIdx, scalarIdx);
                                              });
                       });
}

template<typename NgpSpace>
inline void field_product_impl(const stk::mesh::FieldBase & xField,
                               const stk::mesh::FieldBase & yField,
                               const stk::mesh::FieldBase & zField,
                               const stk::mesh::Selector & selector,
                               const typename NgpSpace::exec_space & execSpace)
{
  field_datatype_execute(xField,
    [&]<typename T>(const stk::mesh::FieldBase&) {
      apply_product_on_field<NgpSpace, T>(xField, yField, zField, selector, execSpace);
    }
  );
}

template<typename NgpSpace, typename Scalar>
inline void field_scale_impl(const Scalar alpha,
                             const stk::mesh::FieldBase & xField,
                             const stk::mesh::Selector & selector,
                             const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto fieldData = xField.data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_scale",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         auto fieldValues = fieldData.bucket_values(bucketIndex);
                         const unsigned numValues = fieldValues.num_entities() * fieldValues.num_scalars();
                         const unsigned numEntities = fieldValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                fieldValues(entityIdx, scalarIdx) *= alpha;
                                              });
                       });
}

template <typename NgpSpace, typename Scalar>
void apply_swap_on_field(const mesh::FieldBase & xField,
                         const mesh::FieldBase & yField,
                         const stk::mesh::Selector & selector,
                         const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto xFieldData = xField.data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace);
  auto yFieldData = yField.data<Scalar, mesh::ReadWrite, NgpSpace>(execSpace);

  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_for("field_swap",
                       TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const TeamHandleType& team) {
                         const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                         auto xValues = xFieldData.bucket_values(bucketIndex);
                         auto yValues = yFieldData.bucket_values(bucketIndex);
                         const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
                         const unsigned numEntities = xValues.num_entities();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numValues),
                                              [&](const unsigned idx) {
                                                const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                Scalar tmp = xValues(entityIdx, scalarIdx);
                                                xValues(entityIdx, scalarIdx) = yValues(entityIdx, scalarIdx);
                                                yValues(entityIdx, scalarIdx) = tmp;
                                              });
                       });
}

template<typename NgpSpace>
inline void field_swap_impl(const stk::mesh::FieldBase & xField,
                            const stk::mesh::FieldBase & yField,
                            const stk::mesh::Selector & selector,
                            const typename NgpSpace::exec_space & execSpace)
{
  field_datatype_execute(xField,
    [&]<typename T>(const stk::mesh::FieldBase&) {
      apply_swap_on_field<NgpSpace, T>(xField, yField, selector, execSpace);
    }
  );
}

template<typename NgpSpace, typename Scalar>
inline void field_dot_impl(Scalar& returnVal,
                           const stk::mesh::FieldBase & xField,
                           const stk::mesh::FieldBase & yField,
                           const stk::mesh::Selector & selector,
                           const typename NgpSpace::exec_space& execSpace)
{
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto xFieldData = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);
  auto yFieldData = yField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);

  auto result = Kokkos::View<Scalar, typename NgpSpace::mem_space>("result");
  auto summer = Kokkos::Sum<Scalar, typename NgpSpace::exec_space>(result);
  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  const unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_reduce("field_dot",
                          TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                          KOKKOS_LAMBDA(const TeamHandleType& team, Scalar& teamReduction) {
                            const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                            auto xValues = xFieldData.bucket_values(bucketIndex);
                            auto yValues = yFieldData.bucket_values(bucketIndex);
                            const unsigned numValues = xValues.num_entities() * xValues.num_scalars();
                            const unsigned numEntities = xValues.num_entities();

                            Scalar localReduction;
                            auto threadSummer = Kokkos::Sum<Scalar, typename NgpSpace::exec_space>(localReduction);
                            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numValues),
                                                    [&](const unsigned idx, Scalar& threadReduction) {
                                                      const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                      const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                      threadReduction += xValues(entityIdx, scalarIdx) * yValues(entityIdx, scalarIdx);
                                                    },
                                                    threadSummer);
                            Kokkos::single(Kokkos::PerTeam(team),
                                           [&]() {
                                             teamReduction += localReduction;
                                           });

                          },
                          summer);
  Kokkos::deep_copy(execSpace, returnVal, result);
  execSpace.fence();

  Scalar globalReturnVal = returnVal;
  stk::all_reduce_sum(xField.get_mesh().parallel(), &returnVal, &globalReturnVal, 1u);
  returnVal = globalReturnVal;
}

template <typename NgpSpace, typename Scalar>
void field_amin_impl(Scalar& aminOut,
                     const stk::mesh::FieldBase & xField,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto fieldData = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);

  auto result = Kokkos::View<Scalar, typename NgpSpace::mem_space>("result");
  auto minner = Kokkos::Min<Scalar, typename NgpSpace::exec_space>(result);
  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  const unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_reduce("field_amin",
                          TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                          KOKKOS_LAMBDA(const TeamHandleType& team, Scalar& teamReduction) {
                            const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                            auto fieldValues = fieldData.bucket_values(bucketIndex);
                            const unsigned numValues = fieldValues.num_entities() * fieldValues.num_scalars();
                            const unsigned numEntities = fieldValues.num_entities();

                            Scalar localReduction;
                            auto threadMinner = Kokkos::Min<Scalar, typename NgpSpace::exec_space>(localReduction);
                            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numValues),
                                                    [&](const unsigned idx, Scalar& threadReduction) {
                                                      const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                      const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                      Scalar fieldValue = fieldValues(entityIdx, scalarIdx);
                                                      Scalar absValue = fieldValue > 0 ? fieldValue : -fieldValue;
                                                      threadReduction = Kokkos::min(threadReduction, absValue);
                                                    },
                                                    threadMinner);
                            Kokkos::single(Kokkos::PerTeam(team),
                                           [&]() {
                                             teamReduction = Kokkos::min(teamReduction, localReduction);
                                           });

                          },
                          minner);
  Kokkos::deep_copy(execSpace, aminOut, result);
  execSpace.fence();

  Scalar localAmin = aminOut;
  stk::all_reduce_min(xField.get_mesh().parallel(), &localAmin, &aminOut, 1u);
}

template <typename NgpSpace, typename Scalar>
void field_asum_impl(Scalar & asumOut,
                     const stk::mesh::FieldBase & xField,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto fieldData = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);

  auto result = Kokkos::View<Scalar, typename NgpSpace::mem_space>("result");
  auto summer = Kokkos::Sum<Scalar, typename NgpSpace::exec_space>(result);
  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  const unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_reduce("field_asum",
                          TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                          KOKKOS_LAMBDA(const TeamHandleType& team, Scalar& teamReduction) {
                            const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                            auto fieldValues = fieldData.bucket_values(bucketIndex);
                            const unsigned numValues = fieldValues.num_entities() * fieldValues.num_scalars();
                            const unsigned numEntities = fieldValues.num_entities();

                            Scalar localReduction;
                            auto threadSummer = Kokkos::Sum<Scalar, typename NgpSpace::exec_space>(localReduction);
                            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numValues),
                                                    [&](const unsigned idx, Scalar& threadReduction) {
                                                      const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                      const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                      Scalar fieldValue = fieldValues(entityIdx, scalarIdx);
                                                      Scalar absValue = fieldValue > 0 ? fieldValue : -fieldValue;
                                                      threadReduction += absValue;
                                                    },
                                                    threadSummer);
                            Kokkos::single(Kokkos::PerTeam(team),
                                           [&]() {
                                             teamReduction += localReduction;
                                           });

                          },
                          summer);
  Kokkos::deep_copy(execSpace, asumOut, result);
  execSpace.fence();

  Scalar localAsum = asumOut;
  stk::all_reduce_sum(xField.get_mesh().parallel(), &localAsum, &asumOut, 1u);
}

template <typename NgpSpace, typename Scalar>
void field_nrm2_impl(Scalar & nrm2Out,
                     const stk::mesh::FieldBase & xField,
                     const stk::mesh::Selector & selector,
                     const typename NgpSpace::exec_space & execSpace)
{
  auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());
  auto fieldData = xField.data<Scalar, mesh::ReadOnly, NgpSpace>(execSpace);

  auto result = Kokkos::View<Scalar, typename NgpSpace::mem_space>("result");
  auto summer = Kokkos::Sum<Scalar, typename NgpSpace::exec_space>(result);
  const stk::topology::rank_t rank = xField.entity_rank();
  const auto& bucketIds = ngpMesh.get_bucket_ids(rank, selector);
  const unsigned numBuckets = bucketIds.size();
  using TeamPolicyType = stk::ngp::TeamPolicy<typename NgpSpace::exec_space>;
  using TeamHandleType = typename TeamPolicyType::member_type;
  Kokkos::parallel_reduce("field_nrm2",
                          TeamPolicyType(execSpace, numBuckets, Kokkos::AUTO),
                          KOKKOS_LAMBDA(const TeamHandleType& team, Scalar& teamReduction) {
                            const int bucketIndex = bucketIds.get<typename NgpSpace::exec_space>(team.league_rank());
                            auto fieldValues = fieldData.bucket_values(bucketIndex);
                            const unsigned numValues = fieldValues.num_entities() * fieldValues.num_scalars();
                            const unsigned numEntities = fieldValues.num_entities();

                            Scalar localReduction;
                            auto threadSummer = Kokkos::Sum<Scalar, typename NgpSpace::exec_space>(localReduction);
                            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numValues),
                                                    [&](const unsigned idx, Scalar& threadReduction) {
                                                      const auto entityIdx = mesh::EntityIdx(idx % numEntities);
                                                      const auto scalarIdx = mesh::ScalarIdx(idx / numEntities);
                                                      Scalar fieldValue = fieldValues(entityIdx, scalarIdx);
                                                      threadReduction += fieldValue * fieldValue;
                                                    },
                                                    threadSummer);
                            Kokkos::single(Kokkos::PerTeam(team),
                                           [&]() {
                                             teamReduction += localReduction;
                                           });

                          },
                          summer);
  Kokkos::deep_copy(execSpace, nrm2Out, result);
  execSpace.fence();

  Scalar localNrm2 = nrm2Out;
  auto comm = xField.get_mesh().parallel();
  stk::all_reduce_sum(comm, &localNrm2, &nrm2Out, 1u);
  nrm2Out = std::sqrt(nrm2Out);
}

//************ end of implementation detail *********************************

} // namespace stk::ngp_field_blas::impl

#endif // STK_MESH_BASEIMPL_NGPFIELDBLASIMPL_HPP
