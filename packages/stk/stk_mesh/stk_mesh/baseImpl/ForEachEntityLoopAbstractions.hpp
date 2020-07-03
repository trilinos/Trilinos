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

#ifndef stk_mesh_ForEachEntityLoopAbstractions_hpp
#define stk_mesh_ForEachEntityLoopAbstractions_hpp

#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, Bucket::size_type, etc
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {
namespace impl {

template <typename ALGORITHM_PER_ENTITY>
void for_each_selected_entity_run(const BulkData &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const ALGORITHM_PER_ENTITY &functor)
{
    const stk::mesh::BucketVector & buckets = mesh.get_buckets(rank, selector);
    const size_t numBuckets = buckets.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t j=0; j<numBuckets; j++)
    {
        stk::mesh::Bucket *bucket = buckets[j];
        for(size_t i=0; i<bucket->size(); i++)
        {
            functor(mesh, stk::mesh::MeshIndex({bucket,i}));
        }
    }
}

template <typename ALGORITHM_PER_ENTITY>
void for_each_selected_entity_run_no_threads(const BulkData &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const ALGORITHM_PER_ENTITY &functor)
{
    const stk::mesh::BucketVector & buckets = mesh.get_buckets(rank, selector);
    for(size_t j=0; j<buckets.size(); j++)
    {
        stk::mesh::Bucket *bucket = buckets[j];
        for(size_t i=0; i<bucket->size(); i++)
        {
            functor(mesh, stk::mesh::MeshIndex({bucket,i}));
        }
    }
}

template <typename ALGORITHM_PER_ENTITY>
void for_each_entity_run_no_threads(const BulkData &mesh, stk::topology::rank_t rank, const ALGORITHM_PER_ENTITY &functor)
{
    const stk::mesh::Selector selectAll = !stk::mesh::Selector();
    for_each_selected_entity_run_no_threads(mesh, rank, selectAll, functor);
}

inline Entity get_entity( const MeshIndex &meshIndex)
{
    return (*meshIndex.bucket)[meshIndex.bucket_ordinal];
}
inline unsigned num_nodes(const MeshIndex &meshIndex)
{
    return meshIndex.bucket->num_nodes(meshIndex.bucket_ordinal);
}
inline Entity const * begin_nodes(const MeshIndex &meshIndex)
{
    return meshIndex.bucket->begin_nodes(meshIndex.bucket_ordinal);
}
inline stk::topology get_topology(const MeshIndex &meshIndex)
{
    return meshIndex.bucket->topology();
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_ForEachEntityLoopAbstractions_hpp

