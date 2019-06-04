// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

#ifndef STK_NGP_REDUCTIONS_H_
#define STK_NGP_REDUCTIONS_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkNgpVector.hpp>

namespace ngp {
template<typename T, typename Index>
KOKKOS_FUNCTION
T reduction_value_type_from_field_value(const T& i, const Index, const T&) 
{
  return i;
}
template<typename T, typename Index>
KOKKOS_FUNCTION
Kokkos::MinMaxScalar<T> reduction_value_type_from_field_value(const T& i, const Index, const Kokkos::MinMaxScalar<T>&)
{
  Kokkos::MinMaxScalar<T> minMax;
  minMax.min_val = i;
  minMax.max_val = i;
  return minMax;
}
template<typename T, typename Index>
KOKKOS_FUNCTION
Kokkos::ValLocScalar<T,Index> reduction_value_type_from_field_value(const T& i, const Index index, const Kokkos::ValLocScalar<T,Index>&)
{
  Kokkos::ValLocScalar<T,Index> val;
  val.val = i;
  val.loc = index;
  return val;
}
template<typename T, typename Index>
KOKKOS_FUNCTION
Kokkos::MinMaxLocScalar<T,Index> reduction_value_type_from_field_value(const T& i, const Index index, const Kokkos::MinMaxLocScalar<T,Index>&)
{
  Kokkos::MinMaxLocScalar<T,Index> val;
  val.min_val = i;
  val.max_val = i;
  val.min_loc = index;
  val.max_loc = index;
  return val;
}
template <typename Mesh, typename Field, typename ReductionOp>
struct ReductionTeamFunctor
{
  using value_type = typename ReductionOp::value_type;
    STK_FUNCTION
    ReductionTeamFunctor(const Mesh m, Field f, stk::NgpVector<unsigned> b, value_type& i)
    : mesh(m), field(f), bucketIds(b), initialValue(i) {}

    using TeamHandleType = typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type;

    STK_FUNCTION
    void operator()(const TeamHandleType& team, value_type& update) const
    {
        const int bucketIndex = bucketIds.device_get(team.league_rank());
        const typename Mesh::BucketType &bucket = mesh.get_bucket(field.get_rank(), bucketIndex);
        unsigned numElements = bucket.size();
        value_type localUpdate = initialValue;
        ReductionOp reduction(localUpdate);
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numElements),
            [&](const int & i, value_type& reduce){
              auto field_value = field.get(typename Mesh::MeshIndex{&bucket, static_cast<unsigned>(i)},0);
              value_type input = reduction_value_type_from_field_value(field_value,i,reduce);
              reduction.join(reduce,input);
            },
            reduction);
        Kokkos::single(Kokkos::PerTeam(team), [&](){reduction.join(update, localUpdate);});
    }

private:
    const Mesh mesh;
    Field field;
    stk::NgpVector<unsigned> bucketIds;
    value_type initialValue;
};

template <typename Mesh, typename Field, typename ReductionOp>
void get_field_reduction(Mesh &mesh, Field field, const stk::mesh::Selector &selector, ReductionOp reduction)
{
    stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(field.get_rank(), selector);
    const unsigned numBuckets = bucketIds.size();
    ReductionTeamFunctor<Mesh, Field, ReductionOp> teamFunctor(mesh, field, bucketIds, reduction.reference());
    Kokkos::parallel_reduce(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO), teamFunctor, reduction);
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


#endif /* STK_NGP_REDUCTIONS_H_ */
