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
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk
{
namespace mesh
{

//void create_interior_block_boundary_sides(BulkData, Selector& blocksToConsider, Part& partToPutSidesInto);

//void create_all_boundary_sides(BulkData, Selector& blocksToSkin, Part& partToPutSidesInto);

void create_exposed_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, Part& partToPutSidesInto)
{
    const PartVector skinnedPart{&partToPutSidesInto};

    ElemElemGraph(bulkData, blocksToSkin).skin_mesh( skinnedPart );
}

Entity get_face_for_element_side_pair(BulkData &bulkData, const SideSetEntry &facet)
{
    const Entity * faces = bulkData.begin_faces(facet.element);
    unsigned numFaces = bulkData.num_faces(facet.element);

    ConnectivityOrdinal const * ordinals = bulkData.begin_face_ordinals(facet.element);

    for(unsigned i = 0; i<numFaces; ++i)
    {
        if(ordinals[i] == facet.side)
            return faces[i];
    }

    return Entity();
}

bool check_global_truth_value(bool truthValue, MPI_Comm communicator)
{
    unsigned localResult = truthValue;
    unsigned globalResult = 0;
    stk::all_reduce_min<unsigned>( communicator, &localResult, &globalResult , 1 );
    return (0 != globalResult);
}

bool check_sideset_inclusion_in_skinning(BulkData &bulkData, stk::mesh::EntityVector &sidesetFaces, Part& skinnedPart)
{
    stk::mesh::EntityVector skinnedFaces;
    stk::mesh::get_selected_entities(skinnedPart, bulkData.buckets(stk::topology::FACE_RANK), skinnedFaces);
    stk::util::sort_and_unique(sidesetFaces);
    return check_global_truth_value(sidesetFaces == skinnedFaces, bulkData.parallel());
}

bool check_if_sideset_is_part_of_skinning(BulkData &bulkData, std::vector<SideSetEntry> &skinnedSideSet, Part& skinnedPart)
{
    for(const SideSetEntry &facet : skinnedSideSet)
    {
        Entity face = get_face_for_element_side_pair(bulkData, facet);

        if(!bulkData.is_valid(face) || !bulkData.bucket(face).member(skinnedPart))
            return false;
    }
    return true;
}

stk::mesh::EntityVector get_faces_from_sideset(BulkData &bulkData, std::vector<SideSetEntry> &skinnedSideSet)
{
    stk::mesh::EntityVector sidesetFaces;

    for(const SideSetEntry &facet : skinnedSideSet)
    {
        Entity face = get_face_for_element_side_pair(bulkData, facet);
        sidesetFaces.push_back(face);
    }

    return sidesetFaces;
}


bool check_exposed_boundary_sides(BulkData &bulkData, const Selector& skinnedBlock, Part& skinnedPart)
{
    std::vector<SideSetEntry> skinnedSideSet = ElemElemGraph(bulkData, skinnedBlock).extract_skinned_sideset(  );

    if(!check_if_sideset_is_part_of_skinning(bulkData, skinnedSideSet, skinnedPart)) return false;

    stk::mesh::EntityVector sidesetFaces = get_faces_from_sideset(bulkData, skinnedSideSet);

    return check_sideset_inclusion_in_skinning(bulkData, sidesetFaces, skinnedPart);
}

} // namespace mesh
} // namespace stk
