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
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>

namespace stk
{
namespace mesh
{

//void create_all_boundary_sides(BulkData &bulkData, const Selector &blocksToConsider, Part &partToPutSidesInto);

//void check_all_boundary_sides(BulkData &bulkData, const Selector &blocksToConsider, Part &partToPutSidesInto);

struct SideSetEntry
{
  SideSetEntry() : element(stk::mesh::Entity()), side(stk::mesh::INVALID_CONNECTIVITY_ORDINAL){};
  SideSetEntry(stk::mesh::Entity in_element, stk::mesh::ConnectivityOrdinal in_side)
    : element(in_element),
      side(in_side)
  {  }

  stk::mesh::Entity element;
  stk::mesh::ConnectivityOrdinal side;
};


class SideSetEntryLess
{
public:
    SideSetEntryLess(const BulkData& mesh);
    bool operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const;
private:
  const BulkData& m_mesh;
};

class SideSetEntryEquals
{
public:
    SideSetEntryEquals(const BulkData& mesh);
    bool operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const;
private:
  const BulkData& m_mesh;
};

//////////////

SideSetEntryLess::SideSetEntryLess(const BulkData& mesh) : m_mesh(mesh){}

bool SideSetEntryLess::operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const
{
    if(m_mesh.identifier(lhs.element) < m_mesh.identifier(rhs.element))
        return true;
    else if(m_mesh.identifier(lhs.element) > m_mesh.identifier(rhs.element))
        return false;
    else
    {
        if(lhs.side<rhs.side)
            return true;
        else
            return false;
    }
    return false;
}

//////////////

SideSetEntryEquals::SideSetEntryEquals(const BulkData& mesh) : m_mesh(mesh){}

bool SideSetEntryEquals::operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const
{
    if(m_mesh.identifier(lhs.element) == m_mesh.identifier(rhs.element) &&
            lhs.side == rhs.side)
        return true;
    return false;
}

class FaceCreator {
public:
    FaceCreator(stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph& elemElemGraph)
    : m_bulkData(bulkData), m_eeGraph(elemElemGraph)
    {

    }

    void fill_side_ordinals(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet, std::vector<int>& ordinals)
    {
        size_t j=element_side_index;
        do
        {
            ordinals.push_back(skinnedSideSet[j].side);
            j++;
        }
        while(j<skinnedSideSet.size() && skinnedSideSet[j].element==skinnedSideSet[element_side_index].element);
    }

    std::vector<int> get_side_ordinals_of_element(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet)
    {
        std::vector<int> ordinals;
        const int max_num_sides_per_any_element = 20;
        ordinals.reserve(max_num_sides_per_any_element);
        fill_side_ordinals(element_side_index, skinnedSideSet, ordinals);
        return ordinals;
    }

    size_t create_face_entities_per_element(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet, const stk::mesh::PartVector& skinParts, std::vector<stk::mesh::sharing_info> &sharedModified)
    {
        std::vector<int> ordinals = get_side_ordinals_of_element(element_side_index, skinnedSideSet);
        stk::mesh::impl::LocalId localId = m_eeGraph.get_local_element_id(skinnedSideSet[element_side_index].element);
        m_eeGraph.create_side_entities(ordinals, localId, skinParts, sharedModified);
        return ordinals.size();
    }

    void create_side_entities_given_sideset(const std::vector<SideSetEntry> &skinnedSideSet, const stk::mesh::PartVector& skinParts, std::vector<stk::mesh::sharing_info>& sharedModified)
    {
        m_bulkData.modification_begin();
        for(size_t i=0;i<skinnedSideSet.size();)
            i += create_face_entities_per_element(i, skinnedSideSet, skinParts, sharedModified);

        std::sort(sharedModified.begin(), sharedModified.end(), SharingInfoLess());
        if(!sharedModified.empty())
        {
            for(size_t i=1; i<sharedModified.size(); i++)
                if(sharedModified[i].m_entity == sharedModified[i-1].m_entity)
                    sharedModified[i].m_owner = sharedModified[i-1].m_owner;
        }

        m_bulkData.make_mesh_parallel_consistent_after_skinning(sharedModified);
    }

private:
    stk::mesh::BulkData& m_bulkData;
    ElemElemGraph& m_eeGraph;
};

class SkinMesh {
public:
    SkinMesh(stk::mesh::BulkData& bulkData, const ElemElemGraph& elemElemGraph, const stk::mesh::PartVector& skinParts, const stk::mesh::Selector& inputSkinSelector, const stk::mesh::Selector* inputAirSelector = nullptr)
    : eeGraph(elemElemGraph), skinSelector(inputSkinSelector), airSelector(inputAirSelector)
    {
    }

    std::vector<int> get_exposed_sides(stk::mesh::impl::LocalId localId, int maxSidesThisElement)
    {
        std::vector<int> exposedSides;
        impl::add_exposed_sides(localId, maxSidesThisElement, eeGraph.get_graph(), exposedSides);
        if(airSelector != nullptr)
            add_exposed_sides_due_to_air_selector(localId, exposedSides);
        return exposedSides;
    }

    void add_exposed_sides_due_to_air_selector(impl::LocalId local_id, std::vector<int> &exposedSides)
    {
        std::vector<int> exposedSidesDueToAirSelector;
        stk::mesh::Entity element = eeGraph.get_entity_from_local_id(local_id);
        std::vector<bool> areAllConnectedElemsAir(eeGraph.get_mesh().bucket(element).topology().num_sides(), true);

        const GraphEdgesForElement& graphEdges = eeGraph.get_edges_for_element(local_id);
        for(const GraphEdge & graphEdge : graphEdges)
        {
            if(is_connected_element_air(graphEdge))
                exposedSidesDueToAirSelector.push_back(graphEdge.side1);
            else
                areAllConnectedElemsAir[graphEdge.side1] = false;
        }

        for(int exposedSide : exposedSidesDueToAirSelector)
            if(areAllConnectedElemsAir[exposedSide])
                exposedSides.push_back(exposedSide);
    }

    bool is_remote_element_air(const ParallelInfoForGraphEdges &parallelInfoForGraphEdges, const stk::mesh::GraphEdge &graphEdge)
    {
        const stk::mesh::impl::parallel_info &parInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
        return parInfo.m_is_air;
    }

    bool is_connected_element_air(const stk::mesh::GraphEdge &graphEdge)
    {
        if(impl::is_local_element(graphEdge.elem2))
            return is_element_selected_and_can_have_side(*airSelector, eeGraph.get_entity_from_local_id(graphEdge.elem2));
        else
            return is_remote_element_air(eeGraph.get_parallel_info_for_graph_edges(), graphEdge);
    }

    bool is_element_selected_and_can_have_side(const stk::mesh::Selector& selector, stk::mesh::Entity otherElement)
    {
        return selector(eeGraph.get_mesh().bucket(otherElement)) && impl::does_element_have_side(eeGraph.get_mesh(), otherElement);
    }

    std::vector<int> get_sides_exposed_on_other_procs(stk::mesh::impl::LocalId localId, int numElemSides)
    {
        std::vector<bool> isConnectedToRemoteElementInBodyToSkin(numElemSides, false);
        std::vector<bool> isOnlyConnectedRemotely(numElemSides, true);
        for(const stk::mesh::GraphEdge &graphEdge : eeGraph.get_edges_for_element(localId))
        {
            if(impl::is_local_element(graphEdge.elem2))
            {
                stk::mesh::Entity other_element = eeGraph.get_entity_from_local_id(graphEdge.elem2);
                if(is_element_selected_and_can_have_side(skinSelector, other_element))
                    isOnlyConnectedRemotely[graphEdge.side1] = false;
            }
            else
            {
                const impl::parallel_info &parallel_edge_info = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
                bool is_other_element_selected = parallel_edge_info.m_in_body_to_be_skinned;
                if(is_other_element_selected)
                    isConnectedToRemoteElementInBodyToSkin[graphEdge.side1] = true;
            }
        }

        std::vector<int> exposedSides;
        for(int side = 0; side < numElemSides; side++)
        {
            if(isConnectedToRemoteElementInBodyToSkin[side] && isOnlyConnectedRemotely[side])
            {
                exposedSides.push_back(side);
            }
        }
        return exposedSides;
    }

    std::vector<int> get_sides_for_skinning(const stk::mesh::Selector& skinSelector,
                                            const stk::mesh::Bucket& bucket,
                                            stk::mesh::Entity element,
                                            stk::mesh::impl::LocalId localId,
                                            const stk::mesh::Selector* airSelector = nullptr)
    {
        int maxSidesThisElement = bucket.topology().num_sides();
        std::vector<int> exposedSides;
        if(stk::mesh::impl::does_element_have_side(bucket.mesh(), element))
        {
            if(skinSelector(bucket))
                exposedSides = get_exposed_sides(localId, maxSidesThisElement);
            else if(airSelector != nullptr && (*airSelector)(bucket))
                exposedSides = get_sides_exposed_on_other_procs(localId, maxSidesThisElement);
        }
        return exposedSides;
    }

    std::vector<SideSetEntry> extract_skinned_sideset()
    {
        std::vector<SideSetEntry> skinnedSideSet;

        const stk::mesh::BulkData& bulk_data = eeGraph.get_mesh();
        const stk::mesh::BucketVector& buckets = bulk_data.get_buckets(stk::topology::ELEM_RANK, bulk_data.mesh_meta_data().locally_owned_part());

        for(size_t i=0;i<buckets.size();++i)
        {
            const stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t j=0;j<bucket.size();++j)
            {
                stk::mesh::Entity element = bucket[j];

                stk::mesh::impl::LocalId localId = eeGraph.get_local_element_id(element);
                std::vector<int> exposedSides = get_sides_for_skinning(skinSelector, bucket, element, localId, airSelector);

                for(size_t k=0; k<exposedSides.size(); ++k)
                {
                    skinnedSideSet.push_back(SideSetEntry(element, static_cast<ConnectivityOrdinal> (exposedSides[k])));
                }
            }
        }

        return skinnedSideSet;
    }

    std::vector<SideSetEntry> extract_interior_sideset()
    {
        std::vector<SideSetEntry> skinnedSideSet;
        const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();
        const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());

        for(const stk::mesh::Bucket* bucket : buckets)
        {
            if(bucket->owned())
            {
                for(size_t i=0;i<bucket->size();++i)
                {
                    stk::mesh::Entity element = (*bucket)[i];
                    impl::LocalId elementId = eeGraph.get_local_element_id(element);
                    ;

                    for(const stk::mesh::GraphEdge & graphEdge : eeGraph.get_edges_for_element(elementId))
                    {
                        stk::mesh::EntityId otherEntityId = 0;
                        stk::mesh::Entity otherElement;
                        bool isElement1InSelector = skinSelector(bucket);
                        bool isElement2InSelector = false;

                        bool isParallelEdge = graphEdge.elem2<0;
                        bool should_add_side = false;

                        if(isParallelEdge)
                        {
                            const impl::parallel_info &parallel_edge_info = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
                            isElement2InSelector = parallel_edge_info.m_in_body_to_be_skinned;
                            if(!isElement1InSelector && !isElement2InSelector) continue;
                            should_add_side = !stk::mesh::impl::are_entity_element_blocks_equivalent(bulkData, element, parallel_edge_info.m_part_ordinals);
                        }
                        else
                        {
                            otherElement = eeGraph.get_entity_from_local_id(graphEdge.elem2);
                            otherEntityId = bulkData.identifier(otherElement);
                            isElement2InSelector = skinSelector(bulkData.bucket(otherElement));
                            if(!isElement1InSelector && !isElement2InSelector) continue;
                            if(!isParallelEdge && bulkData.identifier(element) < otherEntityId)
                                should_add_side = !stk::mesh::impl::are_entity_element_blocks_equivalent(bulkData, element, otherElement);
                        }

                        if(should_add_side)
                        {
                            skinnedSideSet.push_back(SideSetEntry(element, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side1)));
                            if(!isParallelEdge)
                                skinnedSideSet.push_back(SideSetEntry(otherElement, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side2)));
                        }
                    }
                }
            }
        }
        return skinnedSideSet;
    }

private:
    const stk::mesh::ElemElemGraph& eeGraph;
    const stk::mesh::Selector& skinSelector;
    const stk::mesh::Selector* airSelector;
};

///////////////////////////////////////////////////////////////////////////

void create_exposed_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const stk::mesh::PartVector& partToPutSidesInto, stk::mesh::Selector* air)
{
    ElemElemGraph elemElemGraph(bulkData, blocksToSkin, air);
    SkinMesh skinMesh(bulkData, elemElemGraph, partToPutSidesInto, blocksToSkin, air);
    std::vector<SideSetEntry> skinnedSideSet = skinMesh.extract_skinned_sideset();
    stk::util::sort_and_unique(skinnedSideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));
    std::vector<stk::mesh::sharing_info> sharedModified;
    FaceCreator faceCreator(bulkData, elemElemGraph);
    faceCreator.create_side_entities_given_sideset(skinnedSideSet, partToPutSidesInto, sharedModified);
}

void create_exposed_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const stk::mesh::PartVector& partToPutSidesInto)
{
    stk::mesh::Selector *air = nullptr;
    create_exposed_boundary_sides(bulkData, blocksToSkin, partToPutSidesInto, air);
}

void create_interior_block_boundary_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Selector &blocksToConsider, const stk::mesh::PartVector& partToPutSidesInto)
{
    stk::mesh::ElemElemGraph graph(bulkData, blocksToConsider);
    SkinMesh skinMesh(bulkData, graph, partToPutSidesInto, blocksToConsider);
    std::vector<SideSetEntry> skinnedSideSet = skinMesh.extract_interior_sideset();
    stk::util::sort_and_unique(skinnedSideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));
    std::vector<stk::mesh::sharing_info> sharedModified;
    FaceCreator faceCreator(bulkData, graph);
    faceCreator.create_side_entities_given_sideset(skinnedSideSet, partToPutSidesInto, sharedModified);
}

////////////////////////////////////////////////////////////////////////////////////////

Entity get_side_entity_from_ordinal(const std::vector<Entity> &sides, ConnectivityOrdinal const * ordinals, ConnectivityOrdinal requestedOrdinal)
{
    for(unsigned i = 0; i<sides.size(); ++i)
    {
        if(ordinals[i] == requestedOrdinal)
            return sides[i];
    }

    return Entity();
}

Entity get_side_entity_for_element_side_pair(BulkData &bulkData, const SideSetEntry &facet)
{
    const Entity * sides = bulkData.begin(facet.element, bulkData.mesh_meta_data().side_rank());
    ConnectivityOrdinal const * ordinals = bulkData.begin_ordinals(facet.element, bulkData.mesh_meta_data().side_rank());
    std::vector<Entity> sideVector(sides, sides+bulkData.num_sides(facet.element));
    return get_side_entity_from_ordinal(sideVector, ordinals, facet.side);
}

bool check_global_truth_value(bool truthValue, MPI_Comm communicator)
{
    unsigned localResult = truthValue;
    unsigned globalResult = 0;
    stk::all_reduce_min<unsigned>( communicator, &localResult, &globalResult , 1 );
    return (0 != globalResult);
}

stk::mesh::EntityVector get_locally_owned_skinned_sides(BulkData &bulkData, const Part& skinnedPart)
{
    stk::mesh::EntityVector skinnedSides;
    stk::mesh::get_selected_entities(skinnedPart & bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(bulkData.mesh_meta_data().side_rank()), skinnedSides);
    return skinnedSides;
}

bool is_sideset_equivalent_to_skin(BulkData &bulkData, stk::mesh::EntityVector &sidesetSides, const Part& skinnedPart)
{
    stk::mesh::EntityVector skinnedSides = get_locally_owned_skinned_sides(bulkData, skinnedPart);
    stk::util::sort_and_unique(sidesetSides);
    stk::util::sort_and_unique(skinnedSides);
    return check_global_truth_value(sidesetSides == skinnedSides, bulkData.parallel());
}

void add_locally_owned_side_from_element_side_pair(BulkData &bulkData, const SideSetEntry &facet, stk::mesh::EntityVector &sidesetSides)
{
    Entity side = get_side_entity_for_element_side_pair(bulkData, facet);
    if(bulkData.bucket(side).owned())
        sidesetSides.push_back(side);
}

stk::mesh::EntityVector get_locally_owned_sides_from_sideset(BulkData &bulkData, std::vector<SideSetEntry> &skinnedSideSet)
{
    stk::mesh::EntityVector sidesetSides;

    for(const SideSetEntry &facet : skinnedSideSet)
        add_locally_owned_side_from_element_side_pair(bulkData, facet, sidesetSides);

    return sidesetSides;
}


bool check_exposed_boundary_sides(BulkData &bulkData, const Selector& skinnedBlock, Part& skinnedPart)
{
    ElemElemGraph elemElemGraph(bulkData, skinnedBlock);
    SkinMesh skinMesh(bulkData, elemElemGraph, {&skinnedPart}, skinnedBlock);
    std::vector<SideSetEntry> skinnedSideSet = skinMesh.extract_skinned_sideset();

    stk::mesh::EntityVector sidesetSides = get_locally_owned_sides_from_sideset(bulkData, skinnedSideSet);
    return is_sideset_equivalent_to_skin(bulkData, sidesetSides, skinnedPart);
}

bool check_interior_block_boundary_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Selector &skinnedBlock, stk::mesh::Part &skinnedPart)
{
    stk::mesh::ElemElemGraph graph(bulkData, skinnedBlock);
    SkinMesh skinMesh(bulkData, graph, {&skinnedPart}, skinnedBlock);
    std::vector<SideSetEntry> skinnedSideSet = skinMesh.extract_interior_sideset();
    stk::mesh::EntityVector sidesetSides = stk::mesh::get_locally_owned_sides_from_sideset(bulkData, skinnedSideSet);
    return stk::mesh::is_sideset_equivalent_to_skin(bulkData, sidesetSides, skinnedPart);
}


} // namespace mesh
} // namespace stk
