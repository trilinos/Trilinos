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
//

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk
{
namespace mesh
{

//void create_interior_block_boundary_sides(BulkData, Selector& blocksToConsider, Part& partToPutSidesInto);

//void create_all_boundary_sides(BulkData, Selector& blocksToSkin, Part& partToPutSidesInto);

void create_exposed_boundary_sides(BulkData &bulkData, Selector& blocksToSkin, Part& partToPutSidesInto)
{
    const PartVector skinnedPart{&partToPutSidesInto};

    ElemElemGraph elem_elem_graph(bulkData, blocksToSkin);
    elem_elem_graph.skin_mesh( skinnedPart );

}

bool verify_exposed_boundary_sides_count(BulkData &bulkData, const Part& skinnedPart, const std::vector<EntitySidePair> &skinnedSideSet)
{
    std::vector<unsigned> skin_counts;
    count_entities(skinnedPart & bulkData.mesh_meta_data().locally_owned_part(), bulkData, skin_counts);

    std::vector<unsigned>  localSkinnedCount(2);
    localSkinnedCount[0] = skin_counts[bulkData.mesh_meta_data().side_rank()];
    localSkinnedCount[1] = skinnedSideSet.size();
    std::vector<unsigned> globalSkinnedCount = {0, 0};
    stk::all_reduce_sum<unsigned>( bulkData.parallel(), &localSkinnedCount[0], &globalSkinnedCount[0] , 2 );

    return (globalSkinnedCount[0] == globalSkinnedCount[1]);
}

Entity get_face_for_element_side_pair(BulkData &bulkData, const EntitySidePair &facet)
{
    const Entity * faces = bulkData.begin_faces(facet.first);
    unsigned numFaces = bulkData.num_faces(facet.first);

    ConnectivityOrdinal const * ordinals = bulkData.begin_face_ordinals(facet.first);

    for(unsigned i = 0; i<numFaces; ++i)
    {
        if(ordinals[i] == facet.second)
            return faces[i];
    }

    return Entity();
}

bool check_exposed_boundary_sides(BulkData &bulkData, Selector& skinnedBlock, const Part& skinnedPart)
{
    ElemElemGraph elem_elem_graph(bulkData, skinnedBlock);
    std::vector<EntitySidePair> skinnedSideSet = elem_elem_graph.extract_skinned_sideset(  );

    for(const EntitySidePair &facet : skinnedSideSet)
    {
        Entity face = get_face_for_element_side_pair(bulkData, facet);

        if(!bulkData.is_valid(face) || !bulkData.bucket(face).member(skinnedPart))
            return false;
    }

    return verify_exposed_boundary_sides_count(bulkData, skinnedPart, skinnedSideSet);
}


} // namespace mesh
} // namespace stk
