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

#ifndef STK_ELEM_ELEM_GRAPH_IMPL_HPP
#define STK_ELEM_ELEM_GRAPH_IMPL_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include "stk_mesh/base/FEMHelpers.hpp"
#include "GraphTypes.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { struct sharing_info; } }
namespace stk { namespace mesh { class ElemElemGraph; } }
namespace stk { namespace mesh { class Graph; } }
namespace stk { namespace mesh { class ParallelInfoForGraphEdges; } }
namespace stk { namespace mesh { namespace impl { class ElementLocalIdMapper; } } }
namespace stk { class CommSparse; }

namespace stk { namespace mesh {
namespace impl
{

static const LocalId INVALID_LOCAL_ID = std::numeric_limits<impl::LocalId>::max();

struct RemoteElementData
{
public:
    RemoteElementData() : m_other_proc(-1) {}

    RemoteElementData(int proc) :
        m_other_proc(proc) {}

    int get_proc_rank_of_neighbor() const { return m_other_proc; }

    void set_proc_rank(int proc) { m_other_proc = proc; }

    bool operator!=(const RemoteElementData& rhs) const {
        return m_other_proc != rhs.m_other_proc;
    }

private:
    int m_other_proc;
};

struct GraphEdgeProc
{
    GraphEdgeProc(const stk::mesh::EntityId& localElementId, int localSide, const stk::mesh::EntityId &remoteElementId, int remoteSide, int proc_id) :
        m_localElementId(localElementId), m_remoteElementId(remoteElementId), m_localSide(localSide), m_remoteSide(remoteSide), m_proc_id(proc_id) {}
    GraphEdgeProc() : m_localElementId(stk::mesh::InvalidEntityId), m_remoteElementId(stk::mesh::InvalidEntityId), m_localSide(-1),
            m_remoteSide(-1), m_proc_id(-1) {}

    stk::mesh::EntityId get_local_element_global_id() const { return m_localElementId; }
    stk::mesh::EntityId get_remote_element_global_id() const { return m_remoteElementId; }
    int get_local_element_side_index() const { return m_localSide; }
    int get_remote_element_side_index() const { return m_remoteSide; }
    int get_remote_processor_rank() const { return m_proc_id; }

    void set_remote_processor_rank(int proc) { m_proc_id = proc; }

private:
    stk::mesh::EntityId m_localElementId;
    stk::mesh::EntityId m_remoteElementId;
    int m_localSide;
    int m_remoteSide;
    int m_proc_id;
};

struct ParallelInfo
{
public:
    ParallelInfo(int proc, int perm, stk::topology other_elem_topology) :
        m_permutation(perm), m_remote_element_toplogy(other_elem_topology), remoteElementData(proc) {}
    ParallelInfo(int proc, int perm, stk::mesh::EntityId chosen_face_id, stk::topology other_elem_topology) :
        m_permutation(perm), m_remote_element_toplogy(other_elem_topology), remoteElementData(proc) {}

    int get_proc_rank_of_neighbor() const { return remoteElementData.get_proc_rank_of_neighbor(); }

    void set_proc_rank(int proc) { remoteElementData.set_proc_rank(proc); }

    int m_permutation;
    stk::topology m_remote_element_toplogy;

    bool operator!=(const ParallelInfo& rhs) const {
        return m_permutation != rhs.m_permutation ||
                m_remote_element_toplogy != rhs.m_remote_element_toplogy ||
                remoteElementData != rhs.remoteElementData;
    }
private:
    RemoteElementData remoteElementData;
};

inline
std::ostream& operator<<(std::ostream& out, const ParallelInfo& info)
{
    out << "(other_proc=" << info.get_proc_rank_of_neighbor()
            << ", perm=" << info.m_permutation
            << ", remote_top=" << info.m_remote_element_toplogy
            << ")";
    return out;
}

struct SerialElementData
{
public:
    SerialElementData(LocalId elementLocalId, stk::mesh::EntityId elementId, stk::topology elementToplogy, unsigned sideIndex, const stk::mesh::EntityVector& sideNodes) :
        m_elementLocalId(elementLocalId), m_elementIdentifier(elementId), m_elementTopology(elementToplogy), m_sideIndex(sideIndex), m_sideNodes(sideNodes) {}

    SerialElementData()
    : m_elementLocalId(std::numeric_limits<impl::LocalId>::max()),
      m_elementIdentifier(stk::mesh::InvalidEntityId),
      m_elementTopology(stk::topology::INVALID_TOPOLOGY),
      m_sideIndex(std::numeric_limits<unsigned>::max()),
      m_sideNodes(stk::mesh::EntityVector{}) {}

    stk::mesh::EntityId get_element_identifier() const { return m_elementIdentifier; }
    stk::topology get_element_topology() const { return m_elementTopology; }
    const stk::mesh::EntityVector& get_side_nodes() const { return m_sideNodes; }
    LocalId get_element_local_id() const { return m_elementLocalId; }
    unsigned get_element_side_index() const { return m_sideIndex; }
    stk::mesh::Permutation get_permutation() const { return m_perm; }

    bool is_parallel_edge() const { return false; }
    int get_proc_rank_of_neighbor() const { return -1; }

    void clear_side_nodes() { m_sideNodes.clear(); }
    void resize_side_nodes(size_t n) { m_sideNodes.resize(n); }

    void set_element_local_id(LocalId id) { m_elementLocalId = id; }
    void set_element_identifier(stk::mesh::EntityId id) { m_elementIdentifier = id; }
    void set_element_topology(stk::topology topo) { m_elementTopology = topo; }
    void set_element_side_index(unsigned index) { m_sideIndex = index; }
    void set_element_side_nodes(const stk::mesh::EntityVector& nodes) { m_sideNodes = nodes; }
    void set_permutation(stk::mesh::Permutation perm) { m_perm = perm; }

    stk::mesh::EntityVector::iterator side_nodes_begin() { return m_sideNodes.begin(); }

private:
    LocalId m_elementLocalId;
    stk::mesh::EntityId m_elementIdentifier;
    stk::topology m_elementTopology;
    unsigned m_sideIndex;
    stk::mesh::EntityVector m_sideNodes;
    stk::mesh::Permutation m_perm;
};

struct ParallelElementData
{
    ParallelElementData()
    : remoteElementData(),
      serialElementData()
    {}

    bool is_parallel_edge() const { return get_proc_rank_of_neighbor() != -1; }

    stk::mesh::EntityId get_element_identifier() const { return serialElementData.get_element_identifier(); }
    stk::topology get_element_topology() const { return serialElementData.get_element_topology(); }
    const stk::mesh::EntityVector& get_side_nodes() const { return serialElementData.get_side_nodes(); }
    LocalId get_element_local_id() const { return serialElementData.get_element_local_id(); }
    unsigned get_element_side_index() const { return serialElementData.get_element_side_index(); }
    stk::mesh::Permutation get_permutation() const { return serialElementData.get_permutation(); }

    void clear_side_nodes() { serialElementData.clear_side_nodes(); }
    void resize_side_nodes(size_t n) { serialElementData.resize_side_nodes(n); }
    void set_element_local_id(LocalId id) { serialElementData.set_element_local_id(id); }
    void set_element_identifier(stk::mesh::EntityId id) { serialElementData.set_element_identifier(id); }
    void set_element_topology(stk::topology topo) { serialElementData.set_element_topology(topo); }
    void set_element_side_index(unsigned index) { serialElementData.set_element_side_index(index); }
    void set_element_side_nodes(const stk::mesh::EntityVector& nodes) { serialElementData.set_element_side_nodes(nodes); }
    void set_permutation(stk::mesh::Permutation perm) { serialElementData.set_permutation(perm); }

    stk::mesh::EntityVector::iterator side_nodes_begin() { return serialElementData.side_nodes_begin(); }

    int get_proc_rank_of_neighbor() const { return remoteElementData.get_proc_rank_of_neighbor(); }

    void set_proc_rank(int proc) { remoteElementData.set_proc_rank(proc); }

    stk::mesh::EntityId m_suggestedFaceId;

private:
    RemoteElementData remoteElementData;
    SerialElementData serialElementData;
};

struct SharedEdgeInfo
{
public:
    SharedEdgeInfo()
    : m_sharedNodes(),
      m_remoteElementTopology(stk::topology::INVALID_TOPOLOGY),
      remoteElementData(), graphEdgeProc() {}

    stk::mesh::EntityId get_local_element_global_id() const { return graphEdgeProc.get_local_element_global_id(); }
    stk::mesh::EntityId get_remote_element_global_id() const { return graphEdgeProc.get_remote_element_global_id(); }
    int get_local_element_side_index() const { return graphEdgeProc.get_local_element_side_index(); }
    int get_remote_element_side_index() const { return graphEdgeProc.get_remote_element_side_index(); }
    int get_remote_processor_rank() const { return graphEdgeProc.get_remote_processor_rank(); }

    void set_graph_edge_proc(GraphEdgeProc& graph_edge_proc) { graphEdgeProc = graph_edge_proc; }

    int get_proc_rank_of_neighbor() const { return remoteElementData.get_proc_rank_of_neighbor(); }

    void set_proc_rank(int proc) { remoteElementData.set_proc_rank(proc); }

    stk::mesh::EntityVector m_sharedNodes;
    stk::topology m_remoteElementTopology;

private:
    RemoteElementData remoteElementData;
    GraphEdgeProc graphEdgeProc;
};

struct ShellConnectivityData
{
    stk::mesh::EntityId m_nearElementId;
    int                 m_nearElementSide;
    int                 m_nearElementProc;
    stk::mesh::EntityId m_shellElementId;
    stk::mesh::EntityId m_farElementId;
    int                 m_farElementSide;
    int                 m_farElementProc;
    bool                m_farElementIsRemote;
};

struct DeletedElementData
{
    impl::LocalId       m_deletedElement;
    stk::mesh::EntityId m_remoteElement;
    int                 m_remoteProc;
};

struct ElementViaSidePair
{
    stk::mesh::Entity element;
    int side;
};

struct IdViaSidePair
{
    stk::mesh::EntityId id;
    int side;
};

}//namespace impl

const int max_num_sides_per_elem = 10;
const double inverse_of_max_num_sides_per_elem = 0.1;

struct GraphEdge
{
    GraphEdge(impl::LocalId e1, int s1, impl::LocalId e2, int s2)
    {
        set_vertex1(e1,s1);
        set_vertex2(e2,s2);
    }

    GraphEdge() :
        vertex1(std::numeric_limits<impl::LocalId>::max()), vertex2(std::numeric_limits<impl::LocalId>::max())
    {}

    int side1() const { return get_side(vertex1); }
    int side2() const { return get_side(vertex2); }

    void set_vertex1(const impl::LocalId& elem, int side)
    {
        set_vertex(elem, side, vertex1);
    }

    void set_vertex2(const impl::LocalId& elem, int side)
    {
        set_vertex(elem, side, vertex2);
    }

    void set_vertex(const impl::LocalId& elem, int side, impl::LocalId& vertex)
    {
        if(elem>=0)
            vertex = max_num_sides_per_elem*elem+side;
        else
            vertex = max_num_sides_per_elem*elem-side;
    }

    impl::LocalId elem1() const
    {
        return vertex1*inverse_of_max_num_sides_per_elem;
    }

    impl::LocalId elem2() const
    {
        return vertex2*inverse_of_max_num_sides_per_elem;
    }

    int get_side(const impl::LocalId& vertex) const
    {
        return std::abs(vertex)%max_num_sides_per_elem;
    }

    // elem1, side1, elem2, side2

    impl::LocalId vertex1;
    impl::LocalId vertex2;
};

typedef GraphEdge CoincidentElementConnection;

struct GraphEdgeLessByElem2 {
    bool operator()(const GraphEdge& a, const GraphEdge& b) const
    {
        if (a.elem2() != b.elem2())
        {
            return a.elem2() < b.elem2();
        }
        else if (a.side2() != b.side2())
        {
            return a.side2() < b.side2();
        }
        else if (a.elem1() != b.elem1())
        {
            return a.elem1() < b.elem1();
        }
        else
        {
            return a.side1() < b.side1();
        }
    }
};

inline
bool operator==(const GraphEdge& a, const GraphEdge& b)
{
    return  a.vertex1 == b.vertex1 && a.vertex2 == b.vertex2;
}

inline
std::ostream& operator<<(std::ostream& out, const GraphEdge& graphEdge)
{
    out << "(" << graphEdge.vertex1 << " -> " << graphEdge.vertex2 << ")";
    return out;
}

namespace impl {

typedef std::pair<LocalId,int> ElementSidePair;
typedef std::map<GraphEdge, ParallelInfo, GraphEdgeLessByElem2> ParallelGraphInfo;
typedef std::vector<std::vector<LocalId> > ElementGraph;
typedef std::vector<std::vector<int> > SidesForElementGraph;
typedef std::vector<ParallelElementData> ParallelElementDataVector;
typedef std::vector<SerialElementData> SerialElementDataVector;

typedef std::vector<GraphEdge> GraphEdgeVector;

NAMED_PAIR( EntitySidePair , stk::mesh::Entity , entity , unsigned , side_id )
NAMED_PAIR( ProcFaceIdPair , int , proc , stk::mesh::EntityId , side_id )
NAMED_PAIR( ProcVecFaceIdPair , std::vector<int> , proc_vec , stk::mesh::EntityId , side_id )

typedef std::multimap<EntitySidePair, ProcFaceIdPair>  ElemSideToProcAndFaceId;

unsigned get_num_local_elems(const stk::mesh::BulkData& bulkData);

void fill_topologies(const stk::mesh::BulkData& bulkData, const stk::mesh::impl::ElementLocalIdMapper & localMapper, std::vector<stk::topology>& element_topologies);

ElemSideToProcAndFaceId build_element_side_ids_to_proc_map(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &elements_to_communicate);

std::vector<GraphEdgeProc> get_elements_to_communicate(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &killedElements,
        const ElemElemGraph& elem_graph);

std::vector<GraphEdgeProc> communicate_killed_entities(stk::ParallelMachine communicator, const std::vector<GraphEdgeProc>& elements_to_comm);

void pack_elements_to_comm(stk::CommSparse &comm, const std::vector<GraphEdgeProc>& elements_to_comm);

void add_side_into_exposed_boundary(stk::mesh::BulkData& bulkData, const ParallelInfo& parallel_edge_info,
        stk::mesh::Entity local_element, int side_id, stk::mesh::EntityId remote_id, const stk::mesh::PartVector& parts_for_creating_side,
        std::vector<stk::mesh::sharing_info> &shared_modified, stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector, const stk::mesh::PartVector *boundary_mesh_parts = nullptr);

void remove_side_from_death_boundary(stk::mesh::BulkData& bulkData, stk::mesh::Entity local_element,
        stk::mesh::Part &activePart, stk::mesh::EntityVector &deletedEntities, int side_id);

stk::mesh::ConstPartVector get_stk_parts_for_moving_parts_into_death_boundary(const stk::mesh::PartVector *bc_mesh_parts);

int get_element_side_multiplier();

bool is_id_already_in_use_locally(stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank, stk::mesh::EntityId id);

bool does_side_exist_with_different_permutation(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::ConnectivityOrdinal side_ordinal, stk::mesh::Permutation perm);

bool does_element_side_exist(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal);

stk::mesh::Entity connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
        stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts);

void pack_newly_shared_remote_edges(stk::CommSparse &comm, const stk::mesh::BulkData &m_bulk_data, const std::vector<SharedEdgeInfo> &newlySharedEdges);

bool does_element_have_side(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element);

void add_exposed_sides(LocalId elementId, size_t numElemSides, const stk::mesh::Graph& graph, std::vector<int> &element_side_pairs);

void create_sides_created_during_death_part(stk::mesh::MetaData &metaData);
stk::mesh::Part* get_sides_created_during_death_part(const stk::mesh::MetaData &metaData);
void add_parts_from_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::PartVector& side_parts);
stk::mesh::PartVector get_parts_for_creating_side(stk::mesh::BulkData& bulkData, const stk::mesh::PartVector& parts_for_creating_side, stk::mesh::Entity element, int side_ord);
bool side_created_during_death(stk::mesh::BulkData& bulkData, stk::mesh::Entity side);

bool is_local_element(stk::mesh::impl::LocalId elemId);
void fill_element_side_nodes_from_topology(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element, unsigned side_index, stk::mesh::EntityVector& side_nodes);

inline bool is_shell_or_beam2(stk::topology top)
{
    return top.is_shell();
    //return top.is_shell() || top == stk::topology::BEAM_2;
}

}}} // end namespaces stk mesh

#endif
