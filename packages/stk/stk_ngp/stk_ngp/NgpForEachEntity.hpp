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

#ifndef STK_NGP_FOREACHENTITY_H_
#define STK_NGP_FOREACHENTITY_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkNgpVector.hpp>

namespace ngp {

template <typename Mesh, typename AlgorithmPerEntity>
struct ThreadFunctor
{
    STK_FUNCTION
    ThreadFunctor(const typename Mesh::BucketType *b, const AlgorithmPerEntity &f) :
        bucket(b),
        functor(f)
    {}
    STK_FUNCTION
    void operator()(const int& i) const
    {
        functor(typename Mesh::MeshIndex{bucket, static_cast<unsigned>(i)});
    }
    const typename Mesh::BucketType *bucket;
    const AlgorithmPerEntity &functor;
};

template <typename Mesh, typename AlgorithmPerEntity>
struct TeamFunctor
{
    STK_FUNCTION
    TeamFunctor(const Mesh m, const stk::mesh::EntityRank r, stk::NgpVector<unsigned> b, const AlgorithmPerEntity f) :
        mesh(m),
        rank(r),
        bucketIds(b),
        functor(f)
    {}
    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    STK_FUNCTION
    void operator()(const TeamHandleType& team) const
    {
        const int bucketIndex = bucketIds.device_get(team.league_rank());
        const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
        unsigned numElements = bucket.size();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements),
                             ThreadFunctor<Mesh, AlgorithmPerEntity>(&bucket, functor));
    }
    const Mesh mesh;
    const stk::mesh::EntityRank rank;
    stk::NgpVector<unsigned> bucketIds;
    const AlgorithmPerEntity functor;
};

template <typename Mesh, typename AlgorithmPerEntity>
void for_each_entity_run(Mesh &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const AlgorithmPerEntity &functor)
{
    stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(rank, selector);
    unsigned numBuckets = bucketIds.size();
    Kokkos::parallel_for(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO),
                         TeamFunctor<Mesh, AlgorithmPerEntity>(mesh, rank, bucketIds, functor));
}

//    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
//    unsigned numBuckets = mesh.num_buckets(rank);
//    Kokkos::parallel_for(Kokkos::TeamPolicy<MyExecSpace>(numBuckets, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
//    {
//        const int bucketIndex = team.league_rank();
//        const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
//        unsigned numElements = bucket.size();
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&](const int& i)
//        {
//            functor(typename Mesh::MeshIndex{&bucket, static_cast<unsigned>(i)});
//        });
//    });

}


#endif /* STK_NGP_FOREACHENTITY_H_ */
