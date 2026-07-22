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

#ifndef STK_STK_MESH_STK_MESH_BASE_NGPREDUCTIONS
#define STK_STK_MESH_STK_MESH_BASE_NGPREDUCTIONS

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>

namespace stk {
namespace mesh {

template<typename T>
struct identity {
  KOKKOS_FUNCTION
  T operator()(const T t) const {
    return t;
  }
};

template<typename Mesh, typename Field, typename ReductionOp, typename Modifier = identity<typename Field::value_type>>
struct FieldAccessFunctor{
  using value_type = typename ReductionOp::value_type;
  using reduction_op = ReductionOp;
  KOKKOS_FUNCTION
  value_type operator()(const unsigned i, const unsigned j) const
  {
    const int comp_index = (component > -1) ? component : j;
    auto field_value = field.get(typename stk::mesh::FastMeshIndex{bucket_id, static_cast<unsigned>(i)}, comp_index);
    auto value = fm(field_value);
    value_type input = reduction_value_type_from_field_value(value, i, reduction.reference());
    return input;
  }
  KOKKOS_FUNCTION
  stk::mesh::EntityRank get_rank() const { return field.get_rank(); }
  KOKKOS_DEFAULTED_FUNCTION
  FieldAccessFunctor(const FieldAccessFunctor& rhs) = default;
  KOKKOS_FUNCTION
  FieldAccessFunctor(const FieldAccessFunctor& rhs, unsigned b_id, ReductionOp r)
    : field(rhs.field), bucket_id(b_id), reduction(r), fm(Modifier()), component(rhs.component) {}
  KOKKOS_FUNCTION
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

template<typename T>
KOKKOS_FUNCTION
T reduction_value_from_field_value(const T& i, const stk::mesh::EntityId, const T&) 
{
  return i;
}

template<typename T>
KOKKOS_FUNCTION
Kokkos::MinMaxScalar<T> reduction_value_from_field_value(const T& i, const stk::mesh::EntityId, const Kokkos::MinMaxScalar<T>&)
{
  return {i,i};
}

template<typename T>
KOKKOS_FUNCTION
Kokkos::ValLocScalar<T, stk::mesh::EntityId> reduction_value_from_field_value(const T& i, const stk::mesh::EntityId id, const Kokkos::ValLocScalar<T, stk::mesh::EntityId>&)
{
  return {i,id};
}

template<typename T>
KOKKOS_FUNCTION
Kokkos::MinMaxLocScalar<T, stk::mesh::EntityId> reduction_value_from_field_value(const T& i, const stk::mesh::EntityId id, const Kokkos::MinMaxLocScalar<T, stk::mesh::EntityId>&)
{
  return {i,i,id,id};
}

template <typename Mesh, typename AlgorithmPerEntity>
struct ThreadReductionFunctor
{
  KOKKOS_FUNCTION
  ThreadReductionFunctor(const typename Mesh::BucketType *b, const AlgorithmPerEntity &f) :
    bucket(b),
    functor(f)
  {}

  template <typename value_type>
  KOKKOS_FUNCTION
  void operator()(const int& i, value_type & v) const
  {
    functor(typename stk::mesh::FastMeshIndex{bucket->bucket_id(), static_cast<unsigned>(i)}, v);
  }
  const typename Mesh::BucketType *bucket;
  const AlgorithmPerEntity &functor;
};

template <typename Mesh, typename ReductionOpT, typename AlgorithmPerEntity, typename NgpExecSpace>
struct TeamReductionFunctor
{
  using ReductionOp = ReductionOpT;
  using value_type = typename ReductionOp::value_type;
  KOKKOS_FUNCTION
  TeamReductionFunctor(const Mesh m, const stk::mesh::EntityRank r, stk::NgpVector<unsigned> b, ReductionOp & red, const AlgorithmPerEntity& f)
    : mesh(m), rank(r), bucketIds(b), reduction(red), functor(f) {}

  using TeamHandleType = typename stk::ngp::TeamPolicy<NgpExecSpace>::member_type;

  KOKKOS_FUNCTION
  void operator()(const TeamHandleType& team, value_type& team_reduction) const
  {
    const int bucketIndex = bucketIds.get<typename Mesh::MeshExecSpace>(team.league_rank());
    const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
    unsigned numElements = bucket.size();

    value_type my_value;
    reduction.init(my_value);
    ReductionOp thread_reduction(my_value);
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numElements), 
        ThreadReductionFunctor<Mesh, AlgorithmPerEntity>(&bucket, functor), thread_reduction);
    Kokkos::single(Kokkos::PerTeam(team), [&](){
      reduction.join(team_reduction, my_value);
    });
  }
private:
  const Mesh mesh;
  const stk::mesh::EntityRank rank;
  stk::NgpVector<unsigned> bucketIds;
  ReductionOp reduction;
  const AlgorithmPerEntity functor;
};

template <typename Mesh, typename ReductionOp, typename AlgorithmPerEntity>
void for_each_entity_reduce(Mesh& mesh,
    stk::topology::rank_t rank,
    const stk::mesh::Selector& selector,
    ReductionOp& reduction,
    const AlgorithmPerEntity& functor)
{
  stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(rank, selector);
  const unsigned numBuckets = bucketIds.size();
  TeamReductionFunctor<Mesh,ReductionOp, AlgorithmPerEntity, typename Mesh::MeshExecSpace> teamFunctor(mesh, rank, bucketIds, reduction, functor);
  Kokkos::parallel_reduce(stk::ngp::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO), teamFunctor, reduction);
}

template <typename Mesh, typename ReductionOp, typename AlgorithmPerEntity, typename NgpExecSpace>
void for_each_entity_reduce(Mesh& mesh,
    stk::topology::rank_t rank,
    const stk::mesh::Selector& selector,
    ReductionOp& reduction,
    const AlgorithmPerEntity& functor)
{
  stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(rank, selector);
  const unsigned numBuckets = bucketIds.size();
  TeamReductionFunctor<Mesh,ReductionOp, AlgorithmPerEntity, NgpExecSpace> teamFunctor(mesh, rank, bucketIds, reduction, functor);
  Kokkos::parallel_reduce(stk::ngp::TeamPolicy<NgpExecSpace>(numBuckets, Kokkos::AUTO), teamFunctor, reduction);
}

template <typename Mesh, typename Field, typename Reduction, typename Modifier>
struct FieldThreadReductionFunctor
{
  using value_type = typename Reduction::value_type;

  KOKKOS_FUNCTION
  FieldThreadReductionFunctor(Mesh m, Field f, Reduction r, const int c) : mesh(m), field(f), reduction(r), fm(Modifier()), component(c) {}

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex i, value_type & v) const
  {
    const auto min_max = get_min_max_components(i);
      for (int j = min_max.first; j < min_max.second; ++j) {
        auto field_value = fm(field.get(i, j));
        value_type input = reduction_value_from_field_value(field_value, mesh.identifier(mesh.get_entity(field.get_rank(), i)), reduction.reference());
        reduction.join(v, input);
      }
  }

  KOKKOS_FUNCTION
  Kokkos::pair<int,int> get_min_max_components(const stk::mesh::FastMeshIndex i) const
  {
    if (component > -1)
      return Kokkos::make_pair(component, component + 1);
    return Kokkos::make_pair(0, field.get_num_components_per_entity(i));
  }

  const Mesh mesh;
  const Field field;
  const Reduction reduction;
  Modifier fm;
  const int component;
};

//Precondition: Field is synced to device if necessary
template <typename Mesh, typename Field, typename ReductionOp>
void get_field_reduction(Mesh &mesh, Field field, const stk::mesh::Selector &selector, ReductionOp& reduction, const int & component = -1)
{
  for_each_entity_reduce(mesh, field.get_rank(), selector, reduction, FieldThreadReductionFunctor < Mesh, Field,
      ReductionOp, identity<typename Field::value_type>>(mesh, field, reduction, component));
}

//Precondition: Field is synced to device if necessary
template <typename Mesh, typename Field, typename ReductionOp, typename Modifier>
void get_field_reduction(Mesh &mesh, Field field, const stk::mesh::Selector &selector, ReductionOp& reduction, Modifier /*fm*/, const int & component = -1)
{
  for_each_entity_reduce(mesh, field.get_rank(), selector, reduction, FieldThreadReductionFunctor < Mesh, Field,
      ReductionOp, Modifier>(mesh, field, reduction, component));
}

template <typename Mesh, typename Field>
typename Field::value_type get_field_min(Mesh& mesh, Field field, const stk::mesh::Selector& selector)
{
  field.sync_to_device();
  typename Field::value_type reduction_output;
  using Reduction = Kokkos::Min<typename Field::value_type>;
  Reduction min_reduction(reduction_output);
  get_field_reduction(mesh, field, selector, min_reduction);

  return reduction_output;
}

template <typename Mesh, typename Field>
typename Field::value_type get_field_max(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
  field.sync_to_device();
  typename Field::value_type reduction_output;
  using Reduction = Kokkos::Max<typename Field::value_type>;
  Reduction max_reduction(reduction_output);
  get_field_reduction(mesh, field, selector, max_reduction);

  return reduction_output;
}

template <typename Mesh, typename Field>
typename Field::value_type get_field_sum(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
  field.sync_to_device();
  typename Field::value_type reduction_output;
  using Reduction = Kokkos::Sum<typename Field::value_type>;
  Reduction sum_reduction(reduction_output);
  get_field_reduction(mesh, field, selector, sum_reduction);

  return reduction_output;
}

}
}

#endif /* STK_STK_MESH_STK_MESH_BASE_NGPREDUCTIONS */
