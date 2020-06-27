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

#ifndef STK_MESH_NGPREDUCTIONS_HPP
#define STK_MESH_NGPREDUCTIONS_HPP

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/NgpSpaces.hpp>

namespace stk {
namespace mesh {

template<typename T>
STK_FUNCTION
T reduction_value_type_from_field_value(const T& i, const int, const T&) 
{
  return i;
}

template<typename T>
STK_FUNCTION
Kokkos::MinMaxScalar<T> reduction_value_type_from_field_value(const T& i, const unsigned, const Kokkos::MinMaxScalar<T>&)
{
  return {i,i};
}

template<typename T>
STK_FUNCTION
Kokkos::ValLocScalar<T, stk::mesh::EntityId> reduction_value_type_from_field_value(const T& i, const unsigned index, const Kokkos::ValLocScalar<T, stk::mesh::EntityId>&)
{
  return {i,index};
}

template<typename T>
STK_FUNCTION
Kokkos::MinMaxLocScalar<T, stk::mesh::EntityId> reduction_value_type_from_field_value(const T& i, const unsigned index, const Kokkos::MinMaxLocScalar<T, stk::mesh::EntityId>&)
{
  return {i,i,index,index};
}

template<typename T, typename Mesh>
STK_FUNCTION
void replace_loc_with_entity_id(T & t, const Mesh & mesh, const typename Mesh::BucketType & bucket)
{
}

template<typename T, typename Mesh>
STK_FUNCTION
void replace_loc_with_entity_id(Kokkos::ValLocScalar<T, stk::mesh::EntityId> & t, const Mesh & mesh, const typename Mesh::BucketType & bucket)
{
  t.loc = mesh.identifier(bucket[t.loc]);
}

template<typename T, typename Mesh>
STK_FUNCTION
void replace_loc_with_entity_id(Kokkos::MinMaxLocScalar<T, stk::mesh::EntityId> & t, const Mesh & mesh, const typename Mesh::BucketType & bucket)
{
  t.min_loc = mesh.identifier(bucket[t.min_loc]);
  t.max_loc = mesh.identifier(bucket[t.max_loc]);
}

template<typename T>
struct identity {
  STK_FUNCTION
  T operator()(const T t) const {
    return t;
  }
};

template<typename Mesh, typename Field, typename ReductionOp, typename Modifier = identity<typename Field::value_type>>
struct FieldAccessFunctor{
  using value_type = typename ReductionOp::value_type;
  using reduction_op = ReductionOp;
  STK_FUNCTION
  FieldAccessFunctor(Field f, ReductionOp r, const int comp = -1) :
    field(f), bucket_id(0), reduction(r), fm(Modifier()), component(comp) {}
  STK_FUNCTION
  value_type operator()(const unsigned i, const unsigned j) const
  {
    const int comp_index = (component > -1) ? component : j;
    auto field_value = field.get(typename stk::mesh::FastMeshIndex{bucket_id, static_cast<unsigned>(i)}, comp_index);
    auto value = fm(field_value);
    value_type input = reduction_value_type_from_field_value(value, i, reduction.reference());
    return input;
  }
  STK_FUNCTION
  stk::mesh::EntityRank get_rank() const { return field.get_rank(); }
  STK_FUNCTION
  FieldAccessFunctor(const FieldAccessFunctor& rhs) = default;
  STK_FUNCTION
  FieldAccessFunctor(const FieldAccessFunctor& rhs, unsigned b_id, ReductionOp r)
    : field(rhs.field), bucket_id(b_id), reduction(r), fm(Modifier()), component(rhs.component) {}
  STK_FUNCTION
  int num_components(const int i) const {
    if(component > -1) return 1;
    stk::mesh::FastMeshIndex f = {bucket_id, static_cast<unsigned>(i)};
    unsigned nc = field.get_num_components_per_entity(f);
    return nc;
  }
  Field field;
  unsigned bucket_id;
  ReductionOp reduction;
  Modifier fm;
  const int component;
};

template <typename Mesh, typename Accessor>
struct ReductionTeamFunctor
{
  using ReductionOp = typename Accessor::reduction_op;
  using value_type = typename ReductionOp::value_type;
  STK_FUNCTION
  ReductionTeamFunctor(const Mesh m, stk::NgpVector<unsigned> b, Accessor a)
    : mesh(m), bucketIds(b), accessor(a) {}

  using TeamHandleType = typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type;
  STK_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    accessor.reduction.join(dest,src);
  }
  STK_FUNCTION
  void join(volatile value_type& dest, volatile const value_type& src) const {
    accessor.reduction.join(dest,src);
  }

  STK_FUNCTION
  void operator()(const TeamHandleType& team, value_type& update) const
  {
    const int bucketIndex = bucketIds.device_get(team.league_rank());
    const typename Mesh::BucketType &bucket = mesh.get_bucket(accessor.get_rank(), bucketIndex);
    unsigned numElements = bucket.size();
    value_type my_value;
    accessor.reduction.init(my_value);
    ReductionOp reduction(my_value);
    Accessor thread_local_accessor(accessor, bucket.bucket_id(), reduction);
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numElements),
                            [&](unsigned i, value_type& reduce){
      const unsigned nc = thread_local_accessor.num_components(i);
      for(unsigned j = 0; j < nc; ++j){
        value_type input = thread_local_accessor(i,j);
        accessor.reduction.join(reduce,input);
      }
    },
    reduction);
    Kokkos::single(Kokkos::PerTeam(team), [&](){
      replace_loc_with_entity_id(my_value, mesh, bucket);
      accessor.reduction.join(update, my_value);
    });
  }
private:
  const Mesh mesh;
  stk::NgpVector<unsigned> bucketIds;
  Accessor accessor;
};

template <typename Mesh, typename Field, typename ReductionOp>
void get_field_reduction(Mesh &mesh, Field field, const stk::mesh::Selector &selector, ReductionOp& reduction)
{
  stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(field.get_rank(), selector);
  const unsigned numBuckets = bucketIds.size();
  FieldAccessFunctor<Mesh,Field,ReductionOp> accessor(field, reduction);
  ReductionTeamFunctor<Mesh,decltype(accessor)> teamFunctor(mesh, bucketIds, accessor);
  Kokkos::parallel_reduce(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO), teamFunctor, reduction);
}

template <typename Mesh, typename Accessor>
void get_field_reduction(Mesh &mesh, const stk::mesh::Selector &selector, Accessor& accessor)
{
  stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(accessor.get_rank(), selector);
  const unsigned numBuckets = bucketIds.size();
  ReductionTeamFunctor<Mesh,Accessor> teamFunctor(mesh, bucketIds, accessor);
  Kokkos::parallel_reduce(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO), teamFunctor, accessor.reduction);
}

template <typename Mesh, typename Field>
typename Field::value_type get_field_min(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
  field.sync_to_device();
  typename Field::value_type reduction_output;
  Kokkos::Min<typename Field::value_type> min_reduction(reduction_output);
  get_field_reduction(mesh, field, selector, min_reduction);
  return reduction_output;
}

template <typename Mesh, typename Field>
typename Field::value_type get_field_max(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
  field.sync_to_device();
  typename Field::value_type reduction_output;
  Kokkos::Max<typename Field::value_type> max_reduction(reduction_output);
  get_field_reduction(mesh, field, selector, max_reduction);
  return reduction_output;
}

template <typename Mesh, typename Field>
typename Field::value_type get_field_sum(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
  field.sync_to_device();
  typename Field::value_type reduction_output;
  Kokkos::Sum<typename Field::value_type> sum_reduction(reduction_output);
  get_field_reduction(mesh, field, selector, sum_reduction);
  return reduction_output;
}

}
}

#endif
