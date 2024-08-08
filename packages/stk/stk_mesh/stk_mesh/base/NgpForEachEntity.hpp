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

#ifndef STK_MESH_NGP_FOREACHENTITY_HPP
#define STK_MESH_NGP_FOREACHENTITY_HPP

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/DeviceMesh.hpp>
#include <type_traits>

namespace stk {
namespace mesh {

template <typename Mesh, typename AlgorithmPerEntity>
struct ThreadFunctor
{
  KOKKOS_FUNCTION
  ThreadFunctor(const typename Mesh::BucketType *b, const AlgorithmPerEntity &f)
  : bucket(b),
    functor(f)
  {}

  KOKKOS_FUNCTION
  void operator()(const int& i) const
  {
    functor(stk::mesh::FastMeshIndex{bucket->bucket_id(), static_cast<unsigned>(i)});
  }

  const typename Mesh::BucketType *bucket;
  const AlgorithmPerEntity &functor;
};

template <typename Mesh, typename AlgorithmPerEntity, typename EXEC_SPACE>
struct TeamFunctor
{
  using TeamHandleType = typename stk::ngp::TeamPolicy<EXEC_SPACE>::member_type;

  KOKKOS_FUNCTION
  TeamFunctor(const Mesh& m, const stk::mesh::EntityRank r, const stk::NgpVector<unsigned>& b, const AlgorithmPerEntity& f)
  : mesh(m),
    rank(r),
    bucketIds(b),
    functor(f)
  {
  }
 
  KOKKOS_FUNCTION
  void operator()(const TeamHandleType& team) const
  {
    const int bucketIndex = bucketIds.get<EXEC_SPACE>(team.league_rank());
    const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
    unsigned numEntities = bucket.size();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numEntities), ThreadFunctor<Mesh, AlgorithmPerEntity>(&bucket, functor));
  }

  Mesh mesh;
  stk::mesh::EntityRank rank;
  stk::NgpVector<unsigned> bucketIds;
  const AlgorithmPerEntity functor;
};

template <typename Mesh, typename AlgorithmPerEntity>
void for_each_entity_run(Mesh &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const AlgorithmPerEntity &functor)
{
  Kokkos::Profiling::pushRegion("for_each_entity_run with selector");

  stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();

  using EXEC_SPACE = typename Mesh::MeshExecSpace;
  TeamFunctor<Mesh, AlgorithmPerEntity, EXEC_SPACE> teamFunctor(mesh, rank, bucketIds, functor);
  Kokkos::parallel_for(stk::ngp::TeamPolicy<EXEC_SPACE>(numBuckets, Kokkos::AUTO), teamFunctor);

  Kokkos::Profiling::popRegion();
}

template <typename Mesh, typename AlgorithmPerEntity, typename EXEC_SPACE>
void for_each_entity_run(Mesh &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const AlgorithmPerEntity &functor, const EXEC_SPACE& execSpace)
{
  Kokkos::Profiling::pushRegion("for_each_entity_run with selector and EXEC_SPACE");

  stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(rank, selector);
  unsigned numBuckets = bucketIds.size();

  using TeamHandleType = typename stk::ngp::TeamPolicy<EXEC_SPACE>::member_type;
  Kokkos::parallel_for(stk::ngp::TeamPolicy<EXEC_SPACE>(execSpace, numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team){
      const int bucketIndex = bucketIds.get<EXEC_SPACE>(team.league_rank());
      const typename Mesh::BucketType& bucket = mesh.get_bucket(rank, bucketIndex);
      const unsigned numEntities = bucket.size(); 
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numEntities),
        [&](const int& idx) {
          functor(stk::mesh::FastMeshIndex{bucket.bucket_id(), static_cast<unsigned>(idx)});
        }
      );
    }     
  );

  Kokkos::Profiling::popRegion();
}

}
}

#endif
