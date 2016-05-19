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

#include <stk_mesh/base/FEMHelpers.hpp>
#include <sstream>                      // for operator<<, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityId, etc
#include <string>                       // for char_traits, operator<<
#include "stk_mesh/base/MetaData.hpp"   // for get_topology, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsgIf, etc
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

namespace stk
{
namespace mesh
{

namespace
{

void verify_declare_element_side(
        const BulkData & mesh,
        const Entity elem,
        const unsigned local_side_id
        )
{
    stk::topology elem_top = mesh.bucket(elem).topology();

    ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
            "Element[" << mesh.identifier(elem) << "] has no defined topology");

    stk::topology invalid = stk::topology::INVALID_TOPOLOGY;
    stk::topology side_top =
            ((elem_top != stk::topology::INVALID_TOPOLOGY) && (local_side_id < elem_top.num_sides()) )
            ? elem_top.side_topology(local_side_id) : invalid;

    ThrowErrorMsgIf( elem_top!=stk::topology::INVALID_TOPOLOGY && local_side_id >= elem_top.num_sides(),
            "For elem " << mesh.identifier(elem) << ", local_side_id " << local_side_id << ", " <<
            "local_side_id exceeds " << elem_top.name() << ".num_sides() = " << elem_top.num_sides());

    ThrowErrorMsgIf( side_top == stk::topology::INVALID_TOPOLOGY,
            "For elem " << mesh.identifier(elem) << ", local_side_id " << local_side_id << ", " <<
            "No element topology found");
}

void verify_declare_element_edge(
        const BulkData & mesh,
        const Entity elem,
        const unsigned local_edge_id
        )
{
    stk::topology elem_top = mesh.bucket(elem).topology();

    ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
            "Element[" << mesh.identifier(elem) << "] has no defined topology");

    stk::topology invalid = stk::topology::INVALID_TOPOLOGY;
    stk::topology edge_top =
            (elem_top != stk::topology::INVALID_TOPOLOGY && local_edge_id < elem_top.num_edges() )
            ? elem_top.edge_topology() : invalid;

    ThrowErrorMsgIf( elem_top!=stk::topology::INVALID_TOPOLOGY && local_edge_id >= elem_top.num_edges(),
            "For elem " << mesh.identifier(elem) << ", local_edge_id " << local_edge_id << ", " <<
            "local_edge_id exceeds " << elem_top.name() << ".edge_count = " << elem_top.num_edges());

    ThrowErrorMsgIf( edge_top == stk::topology::INVALID_TOPOLOGY,
            "For elem " << mesh.identifier(elem) << ", local_edge_id " << local_edge_id << ", " <<
            "No element topology found");
}

} // unnamed namespace



Entity declare_element(BulkData & mesh,
        PartVector & parts,
        const EntityId elem_id,
        const EntityIdVector & node_ids)
{
    MetaData & fem_meta = MetaData::get(mesh);
    stk::topology top = fem_meta.get_topology(*parts[0]);
    ThrowAssert(node_ids.size() >= top.num_nodes());

    ThrowErrorMsgIf(top == stk::topology::INVALID_TOPOLOGY,
            "Part " << parts[0]->name() << " does not have a local topology");

    PartVector empty;

    const EntityRank entity_rank = stk::topology::ELEMENT_RANK;

    Entity elem = mesh.declare_entity(entity_rank, elem_id, parts);

    const EntityRank node_rank = stk::topology::NODE_RANK;

    Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    PartVector part_scratch;
    part_scratch.reserve(64);

    for(unsigned i = 0; i < top.num_nodes(); ++i)
    {
        //declare node if it doesn't already exist
        Entity node = mesh.get_entity(node_rank, node_ids[i]);
        if(!mesh.is_valid(node))
        {
            node = mesh.declare_entity(node_rank, node_ids[i], empty);
        }

        mesh.declare_relation(elem, node, i, perm, ordinal_scratch, part_scratch);
    }
    return elem;
}

Entity declare_element_to_entity(BulkData & mesh, Entity elem, Entity entity,
        const unsigned relationOrdinal, const PartVector& parts, stk::topology entity_top)
{
    stk::topology elem_top = mesh.bucket(elem).topology();

    std::vector<unsigned> entity_node_ordinals(entity_top.num_nodes());
    elem_top.sub_topology_node_ordinals(mesh.entity_rank(entity), relationOrdinal, entity_node_ordinals.begin());

    const stk::mesh::Entity *elem_nodes = mesh.begin_nodes(elem);
    EntityVector entity_top_nodes(entity_top.num_nodes());
    elem_top.sub_topology_nodes(elem_nodes, mesh.entity_rank(entity), relationOrdinal, entity_top_nodes.begin());

    Permutation perm = mesh.find_permutation(elem_top, elem_nodes, entity_top, &entity_top_nodes[0], relationOrdinal);

    OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    PartVector part_scratch;
    part_scratch.reserve(64);

    stk::mesh::PartVector initialParts;
    initialParts.reserve(parts.size() + 1);
    initialParts = parts;
    initialParts.push_back(&mesh.mesh_meta_data().get_topology_root_part(entity_top));
    mesh.change_entity_parts(entity, initialParts);

    const stk::mesh::ConnectivityOrdinal *side_ordinals = mesh.begin_ordinals(elem, mesh.entity_rank(entity));
    unsigned num_sides = mesh.count_valid_connectivity(elem, mesh.entity_rank(entity));

    bool elem_to_side_exists = false;
    for(unsigned i = 0; i < num_sides; ++i)
    {
        if(side_ordinals[i] == relationOrdinal)
        {
            elem_to_side_exists = true;
            break;
        }
    }

    if(!elem_to_side_exists)
    {
        mesh.declare_relation(elem, entity, relationOrdinal, perm, ordinal_scratch, part_scratch);
    }

    const unsigned num_side_nodes = mesh.count_valid_connectivity(entity, stk::topology::NODE_RANK);
    if(0 == num_side_nodes)
    {
        Permutation node_perm = stk::mesh::Permutation::INVALID_PERMUTATION;
        Entity const *elem_nodes_local = mesh.begin_nodes(elem);
        for(unsigned i = 0; i < entity_top.num_nodes(); ++i)
        {
            Entity node = elem_nodes_local[entity_node_ordinals[i]];
            mesh.declare_relation(entity, node, i, node_perm, ordinal_scratch, part_scratch);
        }
    }
    else
    {
        ThrowAssertMsg(num_side_nodes == entity_top.num_nodes(),
                "declare_element_to_entity: " << mesh.entity_key(entity) << " already exists with different number of nodes.");
    }

    return entity;
}


Entity declare_element_side(
        BulkData & mesh,
        Entity elem,
        Entity side,
        const unsigned local_side_id,
        const stk::mesh::PartVector& parts)
{
    verify_declare_element_side(mesh, elem, local_side_id);

    stk::topology elem_top = mesh.bucket(elem).topology();
    stk::topology side_top = elem_top.side_topology(local_side_id);
    return declare_element_to_entity(mesh, elem, side, local_side_id, parts, side_top);
}

stk::mesh::EntityId side_id_formula(const stk::mesh::BulkData& bulkData, stk::mesh::Entity elem, unsigned sideOrdinal)
{
    //this is the side-id formula used by IO. the "+1" is because IO always uses one-based side ordinals
    return 10*bulkData.identifier(elem) + sideOrdinal+1;
}

Entity declare_element_side(BulkData& bulkData, Entity elem, const unsigned side_ordinal, const stk::mesh::PartVector& add_parts)
{
    stk::mesh::Entity sideEntity = stk::mesh::impl::get_side_for_element(bulkData, elem, side_ordinal);
    if(bulkData.is_valid(sideEntity))
    {
        bulkData.change_entity_parts(sideEntity, add_parts, {});
    }
    else
    {
        stk::mesh::EntityId global_side_id = side_id_formula(bulkData, elem, side_ordinal);
        sideEntity = stk::mesh::declare_element_side(bulkData, global_side_id, elem, side_ordinal, add_parts);
        if (bulkData.has_face_adjacent_element_graph()) {
            stk::mesh::ElemElemGraph &graph = bulkData.get_face_adjacent_element_graph();
            stk::mesh::SideConnector sideConnector = graph.get_side_connector();
            sideConnector.connect_side_to_all_elements(sideEntity, elem, side_ordinal);
        }
    }
    return sideEntity;
}

Entity declare_element_side( BulkData & mesh ,
                             const stk::mesh::EntityId global_side_id ,
                             Entity elem ,
                             const unsigned local_side_id ,
                             stk::mesh::Part* part)
{
    stk::mesh::PartVector parts;
    if(part!=NULL)
    {
        parts.push_back(part);
    }
    return declare_element_side(mesh, global_side_id, elem, local_side_id, parts);
}

Entity declare_element_side( BulkData & mesh ,
                               Entity elem ,
                               Entity side ,
                               const unsigned local_side_id ,
                               stk::mesh::Part* part)
{
    stk::mesh::PartVector parts;
    if(part!=NULL)
    {
        parts.push_back(part);
    }
    return declare_element_side(mesh, elem, side, local_side_id, parts);
}

Entity declare_element_edge(
        BulkData & mesh,
        Entity elem,
        Entity edge,
        const unsigned local_edge_id,
        const stk::mesh::PartVector& parts)
{
    verify_declare_element_edge(mesh, elem, local_edge_id);
    stk::topology elem_top = mesh.bucket(elem).topology();
    stk::topology edge_top = elem_top.edge_topology();
    return declare_element_to_entity(mesh, elem, edge, local_edge_id, parts, edge_top);
}

Entity declare_element_side(
        BulkData & mesh,
        const stk::mesh::EntityId global_side_id,
        Entity elem,
        const unsigned local_side_id,
        const stk::mesh::PartVector& parts)
{
    verify_declare_element_side(mesh, elem, local_side_id);

    stk::topology elem_top = mesh.bucket(elem).topology();
    stk::topology side_top = elem_top.side_topology(local_side_id);

    PartVector empty_parts;
    Entity side = mesh.get_entity(side_top.rank(), global_side_id);
    if(!mesh.is_valid(side))
    {
        side = mesh.declare_entity(side_top.rank(), global_side_id, empty_parts);
// It seem like declare_element_side should be called even if the side is valid to make sure it is attached to the element.
        declare_element_side(mesh, elem, side, local_side_id, parts);
    }
    return side;
}

Entity declare_element_edge(
        BulkData & mesh,
        const stk::mesh::EntityId global_edge_id,
        Entity elem,
        const unsigned local_edge_id,
        const stk::mesh::PartVector &parts)
{
    verify_declare_element_edge(mesh, elem, local_edge_id);

    stk::topology elem_top = mesh.bucket(elem).topology();
    stk::topology edge_top = elem_top.edge_topology();

    PartVector empty_parts;
    Entity edge = mesh.get_entity(edge_top.rank(), global_edge_id);
    if(!mesh.is_valid(edge))
    {
        edge = mesh.declare_entity(edge_top.rank(), global_edge_id, empty_parts);
// It seem like declare_element_edge should be called even if the edge is valid to make sure it is attached to the element.
        declare_element_edge(mesh, elem, edge, local_edge_id, parts);
    }
    return edge;
}

typedef std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ConnectivityAndOrdinal;


void mark_if_negative_permutation(std::pair<bool, unsigned> &result, stk::topology sub_topology)
{
    if (result.first && result.second >= sub_topology.num_positive_permutations())
    {
        result.first = false;
    }
}

void set_ordinal_and_permutation_if_equivalent(std::pair<bool, unsigned> &result,
                                               unsigned ordinal,
                                               ConnectivityAndOrdinal &ordinalAndPermutation)
{
    if (result.first == true)
    {
        ordinalAndPermutation.first = static_cast<stk::mesh::ConnectivityOrdinal>(ordinal);
        ordinalAndPermutation.second = static_cast<stk::mesh::Permutation>(result.second);
    }
}


class ShellPermutationFilter
{
public:
    ShellPermutationFilter(const stk::mesh::BulkData& mesh,
                           stk::mesh::Entity parent_entity,
                           stk::mesh::EntityRank to_rank)
    : m_mesh(mesh), m_entity(parent_entity), m_toRank(to_rank), m_filterForShell(false)
    {
        stk::topology elemTopology = m_mesh.bucket(m_entity).topology();
        m_filterForShell = elemTopology.is_shell() && to_rank == mesh.mesh_meta_data().side_rank();
    }

    ~ShellPermutationFilter() {}

    bool set_ordinal_and_permutation(const stk::mesh::EntityVector &nodes_of_sub_rank,
                                     stk::mesh::Entity nodes_of_sub_topology[],
                                     stk::topology sub_topology,
                                     unsigned ordinal,
                                     ConnectivityAndOrdinal &ordinalAndPermutation) const
    {
        std::pair<bool, unsigned> result = sub_topology.equivalent(nodes_of_sub_rank, nodes_of_sub_topology);

        if(m_filterForShell)
            mark_if_negative_permutation(result, sub_topology);

        set_ordinal_and_permutation_if_equivalent(result, ordinal, ordinalAndPermutation);
        return result.first;
    }

private:
    ShellPermutationFilter();

    const stk::mesh::BulkData& m_mesh;
    stk::mesh::Entity m_entity;
    stk::mesh::EntityRank m_toRank;
    bool m_filterForShell;
};

class ActivePermutationFilter
{
public:
    ActivePermutationFilter() {}
    ~ActivePermutationFilter() {}

    bool set_ordinal_and_permutation(const stk::mesh::EntityVector &nodes_of_sub_rank,
                                     stk::mesh::Entity nodes_of_sub_topology[],
                                     stk::topology sub_topology,
                                     unsigned ordinal,
                                     ConnectivityAndOrdinal &ordinalAndPermutation) const
    {
        std::pair<bool, unsigned> result = sub_topology.equivalent(nodes_of_sub_rank, nodes_of_sub_topology);
        mark_if_negative_permutation(result, sub_topology);
        set_ordinal_and_permutation_if_equivalent(result, ordinal, ordinalAndPermutation);
        return result.first;
    }
};

template<typename PermutationFilter>
OrdinalAndPermutation
get_ordinal_and_permutation_with_filter(const stk::mesh::BulkData& mesh,
                                        stk::mesh::Entity parent_entity,
                                        stk::mesh::EntityRank to_rank,
                                        const stk::mesh::EntityVector &nodes_of_sub_rank,
                                        PermutationFilter &pFilter)
{
    ConnectivityAndOrdinal ordinalAndPermutation = std::make_pair(stk::mesh::INVALID_CONNECTIVITY_ORDINAL, stk::mesh::INVALID_PERMUTATION);

    unsigned nodes_of_sub_rank_size = nodes_of_sub_rank.size();

    const Entity* elemNodes = mesh.begin_nodes(parent_entity);
    stk::topology elemTopology = mesh.bucket(parent_entity).topology();
    unsigned num_entities_of_sub_topology = elemTopology.num_sub_topology(to_rank);
    const unsigned max_nodes_possible = 100;
    stk::mesh::Entity nodes_of_sub_topology[max_nodes_possible];

    for (unsigned i=0;i<num_entities_of_sub_topology;++i)
    {
        stk::topology sub_topology = elemTopology.sub_topology(to_rank, i);
        unsigned num_nodes = sub_topology.num_nodes();

        if (num_nodes !=  nodes_of_sub_rank_size)
        {
            continue;
        }

        ThrowRequireMsg(num_nodes<=max_nodes_possible, "Program error. Exceeded expected array dimensions. Contact sierra-help for support.");
        elemTopology.sub_topology_nodes(elemNodes, to_rank, i, nodes_of_sub_topology);

        pFilter.set_ordinal_and_permutation(nodes_of_sub_rank, nodes_of_sub_topology, sub_topology, i, ordinalAndPermutation);
        if (ordinalAndPermutation.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL) {
            break;
        }
    }

    return ordinalAndPermutation;
}

OrdinalAndPermutation
get_ordinal_and_permutation(const stk::mesh::BulkData& mesh,
                            stk::mesh::Entity parent_entity,
                            stk::mesh::EntityRank to_rank,
                            const stk::mesh::EntityVector &nodes_of_sub_rank)
{
    ShellPermutationFilter pFilter(mesh, parent_entity, to_rank);
    return get_ordinal_and_permutation_with_filter(mesh, parent_entity, to_rank, nodes_of_sub_rank, pFilter);
}

OrdinalAndPermutation
get_ordinal_and_positive_permutation(const stk::mesh::BulkData& mesh,
                                     stk::mesh::Entity parent_entity,
                                     stk::mesh::EntityRank to_rank,
                                     const stk::mesh::EntityVector &nodes_of_sub_rank)
{
    ActivePermutationFilter pFilter;
    return get_ordinal_and_permutation_with_filter(mesh, parent_entity, to_rank, nodes_of_sub_rank, pFilter);
}

stk::topology get_subcell_nodes(const BulkData& mesh, const Entity entity,
        EntityRank subcell_rank,
        unsigned subcell_identifier,
        EntityVector & subcell_nodes)
{
    ThrowAssert(subcell_rank <= stk::topology::ELEMENT_RANK);

    subcell_nodes.clear();

    // get cell topology
    stk::topology celltopology = mesh.bucket(entity).topology();

    //error checking
    {
//no celltopology defined
        if(celltopology == stk::topology::INVALID_TOPOLOGY)
        {
            return celltopology;
        }

// valid ranks fall within the dimension of the cell topology
        const bool bad_rank = subcell_rank >= celltopology.dimension();
        ThrowInvalidArgMsgIf( bad_rank, "subcell_rank is >= celltopology dimension\n");

// subcell_identifier must be less than the subcell count
        const bool bad_id = subcell_identifier >= celltopology.num_sub_topology(subcell_rank);
        ThrowInvalidArgMsgIf( bad_id, "subcell_id is >= subcell_count\n");
    }

    // Get the cell topology of the subcell
    stk::topology subcell_topology =
            celltopology.sub_topology(subcell_rank, subcell_identifier);

    const int num_nodes_in_subcell = subcell_topology.num_nodes();

    // For the subcell, get it's local nodes ids
    std::vector<unsigned> subcell_node_local_ids(num_nodes_in_subcell);
    celltopology.sub_topology_node_ordinals(subcell_rank, subcell_identifier, subcell_node_local_ids.begin());

    Entity const *node_relations = mesh.begin_nodes(entity);
    subcell_nodes.reserve(num_nodes_in_subcell);

    for(int i = 0; i < num_nodes_in_subcell; ++i)
    {
        subcell_nodes.push_back(node_relations[subcell_node_local_ids[i]]);
    }

    return subcell_topology;
}

int get_entity_subcell_id(const BulkData& mesh,
        const Entity entity,
        const EntityRank subcell_rank,
        stk::topology subcell_topology,
        const std::vector<Entity>& subcell_nodes)
{
    ThrowAssert(subcell_rank <= stk::topology::ELEMENT_RANK);
    const int INVALID_SIDE = -1;

    stk::topology entity_topology = mesh.bucket(entity).topology();
    ThrowAssert(entity_topology.num_nodes() == mesh.num_nodes(entity));

    const Entity* entity_nodes = mesh.begin_nodes(entity);
    std::vector<Entity> topology_sub_nodes(subcell_nodes.size());

    for(size_t i = 0; i < entity_topology.num_sub_topology(subcell_rank); ++i)
    {
        entity_topology.sub_topology_nodes(entity_nodes, subcell_rank, i, topology_sub_nodes.begin());
        if(subcell_topology.equivalent(topology_sub_nodes, subcell_nodes).first)
        {
            return i;
        }
    }

    return INVALID_SIDE;
}

}
}
