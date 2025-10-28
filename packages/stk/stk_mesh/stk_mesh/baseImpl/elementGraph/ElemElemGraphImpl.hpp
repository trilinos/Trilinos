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

#ifndef STK_ELEM_ELEM_GRAPH_IMPL_HPP
#define STK_ELEM_ELEM_GRAPH_IMPL_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/SortAndUnique.hpp>
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

struct RemoteElementData
{
public:
    RemoteElementData() : m_other_proc(-1) {}

    RemoteElementData(int proc) :
        m_other_proc(proc) {}

    int get_proc_rank_of_neighbor() const { return m_other_proc; }

    void set_proc_rank(int proc) { m_other_proc = proc; }

    bool operator==(const RemoteElementData& rhs) const {
        return m_other_proc == rhs.m_other_proc;
    }

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
        m_permutation(perm), m_remote_element_topology(other_elem_topology), remoteElementData(proc) {}
    ParallelInfo(int proc, int perm, stk::mesh::EntityId /*chosen_face_id*/, stk::topology other_elem_topology) :
        m_permutation(perm), m_remote_element_topology(other_elem_topology), remoteElementData(proc) {}
    ParallelInfo() :
        m_permutation(INVALID_PERMUTATION), m_remote_element_topology(stk::topology::INVALID_TOPOLOGY), remoteElementData(-1) {}

    int get_proc_rank_of_neighbor() const { return remoteElementData.get_proc_rank_of_neighbor(); }

    void set_proc_rank(int proc) { remoteElementData.set_proc_rank(proc); }

    int m_permutation;
    stk::topology m_remote_element_topology;

    bool operator==(const ParallelInfo& rhs) const {
        return m_permutation == rhs.m_permutation && 
                m_remote_element_topology == rhs.m_remote_element_topology &&
                remoteElementData == rhs.remoteElementData;
    }

    bool operator!=(const ParallelInfo& rhs) const {
        return m_permutation != rhs.m_permutation ||
                m_remote_element_topology != rhs.m_remote_element_topology ||
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
            << ", remote_top=" << info.m_remote_element_topology
            << ")";
    return out;
}

struct SerialElementData
{
public:
    SerialElementData(LocalId elementLocalId, stk::mesh::EntityId elementId, stk::topology elementTopology, unsigned sideIndex, const stk::mesh::EntityVector& sideNodes) :
        m_elementLocalId(elementLocalId), m_elementIdentifier(elementId), m_elementTopology(elementTopology), m_sideIndex(sideIndex), m_sideNodes(sideNodes) {}

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

    stk::mesh::Entity * side_nodes_begin() { return m_sideNodes.data(); }

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

    stk::mesh::Entity * side_nodes_begin() { return serialElementData.side_nodes_begin(); }

    int get_proc_rank_of_neighbor() const { return remoteElementData.get_proc_rank_of_neighbor(); }

    void set_proc_rank(int proc) { remoteElementData.set_proc_rank(proc); }

    stk::mesh::EntityId m_suggestedFaceId{0u};

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

struct ConnectedElementInfo
{
    stk::mesh::Entity element;
    int thisSide;
    int otherSide;
};

struct IdViaSidePair
{
    stk::mesh::EntityId id;
    int side;
};

}//namespace impl

constexpr int max_num_sides_per_elem = 8;

struct GraphEdge
{
    GraphEdge(impl::LocalId e1, int s1, impl::LocalId e2, int s2)
    {
        set_vertex1(e1,s1);
        set_vertex2(e2,s2);
    }

    GraphEdge() :
        vertex1(impl::INVALID_LOCAL_ID), vertex2(impl::INVALID_LOCAL_ID)
    {}

    GraphEdge(const GraphEdge& rhs)
    : vertex1(rhs.vertex1), vertex2(rhs.vertex2)
    {}

    GraphEdge(const GraphEdge&& rhs)
    : vertex1(rhs.vertex1), vertex2(rhs.vertex2)
    {}

    GraphEdge& operator=(const GraphEdge&) = default;
    GraphEdge& operator=(GraphEdge&&) = default;

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
        return vertex1/max_num_sides_per_elem;
    }

    impl::LocalId elem2() const
    {
        return vertex2/max_num_sides_per_elem;
    }

    int get_side(const impl::LocalId& vertex) const
    {
        return std::abs(vertex)%max_num_sides_per_elem;
    }

    bool is_elem2_local() const
    {
      return vertex2 >= 0;
    }

    // elem1, side1, elem2, side2

    impl::LocalId vertex1;
    impl::LocalId vertex2;
};

constexpr bool is_valid(const GraphEdge& lhs)
{
    return lhs.vertex1 != impl::INVALID_LOCAL_ID;
}

using CoincidentElementConnection = GraphEdge;

struct GraphEdgeLessByElem1 {
    bool operator()(const GraphEdge& a, const GraphEdge& b) const
    {
        impl::LocalId a_elem1 = a.elem1();
        impl::LocalId b_elem1 = b.elem1();

        if (a_elem1 != b_elem1)
        {
            return a_elem1 < b_elem1;
        }

        impl::LocalId a_elem2 = std::abs(a.elem2());
        impl::LocalId b_elem2 = std::abs(b.elem2());
        if (a_elem2 != b_elem2)
        {
            return a_elem2 < b_elem2;
        }

        int a_side2 = a.side2();
        int b_side2 = b.side2();
        if (a_side2 != b_side2)
        {
            return a_side2 < b_side2;
        }

        int a_side1 = a.side1();
        int b_side1 = b.side1();
        if(a_side1 == b_side1)
        {
          return a.elem2() < b.elem2();
        }

        return a_side1 < b_side1;
    }
};

struct GraphEdgeLessByElem2Only
{
    bool operator()(const GraphEdge& a, const GraphEdge& b) const
    {
        impl::LocalId a_elem2 = std::abs(a.elem2());
        impl::LocalId b_elem2 = std::abs(b.elem2());

        return a_elem2 < b_elem2 || (a_elem2 == b_elem2 && a.side2() < b.side2());
    }  
};

inline
bool operator<(const GraphEdge& a, const GraphEdge& b)
{
  GraphEdgeLessByElem1 lessByElem1;
  return lessByElem1(a, b);
}

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

    bool operator()(const std::pair<GraphEdge,impl::ParallelInfo>& a, const GraphEdge& b) const
    {
        return operator()(a.first, b);
    }
    bool operator()(const GraphEdge& a, const std::pair<GraphEdge,impl::ParallelInfo>& b) const
    {
        return operator()(a, b.first);
    }
    bool operator()(const std::pair<GraphEdge,impl::ParallelInfo>& a, const std::pair<GraphEdge,impl::ParallelInfo>& b) const
    {
        return operator()(a.first, b.first);
    }
};

inline
bool operator==(const GraphEdge& a, const GraphEdge& b)
{
    return  a.vertex1 == b.vertex1 && a.vertex2 == b.vertex2;
}

inline
bool operator!=(const GraphEdge& a, const GraphEdge& b)
{
    return  !(a == b);
}

inline
std::ostream& operator<<(std::ostream& out, const GraphEdge& graphEdge)
{
    out << "GraphEdge vertices: (" << graphEdge.vertex1 << " -> " << graphEdge.vertex2 
        << "), element-side pairs: (" << graphEdge.elem1() << ", " << graphEdge.side1() 
        << ") -> (" << graphEdge.elem2() << ", " << graphEdge.side2() << ")";
    return out;
}

namespace impl {

typedef std::pair<LocalId,int> ElementSidePair;
typedef std::vector<std::pair<GraphEdge,ParallelInfo>> ParallelGraphInfo;
typedef std::vector<std::vector<LocalId> > ElementGraph;
typedef std::vector<std::vector<int> > SidesForElementGraph;
typedef std::vector<ParallelElementData> ParallelElementDataVector;
typedef std::vector<SerialElementData> SerialElementDataVector;

typedef std::vector<GraphEdge> GraphEdgeVector;

NAMED_PAIR( EntitySidePair , stk::mesh::Entity , entity , unsigned , side_id )

struct ElemSideProc{
  ElemSideProc(stk::mesh::Entity elem, unsigned side, int procArg)
  : elemSidePair(elem, side), proc(procArg)
  {}

  EntitySidePair elemSidePair;
  int proc;

  bool operator<(const EntitySidePair& rhs) const
  {
    return elemSidePair < rhs;
  }

  bool operator<(const ElemSideProc& rhs) const
  {
    return elemSidePair < rhs.elemSidePair;
  }
};

inline
bool operator<(const EntitySidePair& lhs, const ElemSideProc& rhs)
{
  return lhs < rhs.elemSidePair;
}

using ElemSideProcVector = std::vector<ElemSideProc>;

unsigned get_num_local_elems(const stk::mesh::BulkData& bulkData);

bool fill_topologies(stk::mesh::ElemElemGraph& eeGraph, const stk::mesh::impl::ElementLocalIdMapper & localMapper, std::vector<stk::topology>& element_topologies);

ElemSideProcVector build_element_side_ids_to_proc_map(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &elements_to_communicate);

std::vector<GraphEdgeProc> get_elements_to_communicate(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &killedElements,
        const ElemElemGraph& elem_graph);

std::vector<GraphEdgeProc> communicate_killed_entities(stk::ParallelMachine communicator, const std::vector<GraphEdgeProc>& elements_to_comm);

void pack_elements_to_comm(stk::CommSparse &comm, const std::vector<GraphEdgeProc>& elements_to_comm);

bool can_add_side_into_exposed_boundary(const stk::mesh::BulkData& bulkData, stk::mesh::Entity local_element, int side_id,
                                        const stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                                        const stk::mesh::Part& activePart);

void add_side_into_exposed_boundary(stk::mesh::BulkData& bulkData, const ParallelInfo& parallel_edge_info,
        stk::mesh::Entity local_element, int side_id, const stk::mesh::PartVector& parts_for_creating_side,
        std::vector<stk::mesh::sharing_info> &shared_modified, stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
        const stk::mesh::Part& activePart, const stk::mesh::PartVector *boundary_mesh_parts = nullptr);

int get_number_of_connected_active_elements(const stk::mesh::BulkData& bulkData,
                                            stk::mesh::Entity localElement,
                                            int localOrdinal,
                                            const stk::mesh::Part& activePart,
                                            const stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector);

bool is_exposed_side(const stk::mesh::BulkData& bulkData, stk::mesh::Entity local_element, int side_id,
                     const stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                     const stk::mesh::Part& activePart);

bool remove_side_from_death_boundary(stk::mesh::BulkData& bulkData, stk::mesh::Entity local_element,
        stk::mesh::Part &activePart, stk::mesh::EntityVector &deletedEntities, int side_id, stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector);

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
void fill_element_side_nodes_from_topology(stk::topology elemTopo, const stk::mesh::Entity* elemNodes, unsigned side_index, stk::mesh::EntityVector& side_nodes);

inline bool is_shell_or_beam2(stk::topology top)
{
    return top.is_shell();
    //return top.is_shell() || top == stk::topology::BEAM_2;
}

inline bool is_shell_face_side(stk::topology elemTopology, unsigned sideIndex)
{
  return elemTopology.is_shell() && (sideIndex == 0 || elemTopology.side_topology(sideIndex) == elemTopology.side_topology(0));
}

struct TopologyChecker
{
    bool are_both_shells() const
    {
        return is_shell_or_beam2(localTopology) && is_shell_or_beam2(remoteTopology);
    }

    bool are_both_not_shells() const
    {
        return !is_shell_or_beam2(localTopology) && !is_shell_or_beam2(remoteTopology);
    }

    stk::topology localTopology;
    stk::topology remoteTopology;
};

template<typename Key, typename Value>
struct ComparePairByKey {

  inline bool operator()(const std::pair<Key, Value>& lhs, const std::pair<Key, Value>& rhs) const
  {
    return lhs.first < rhs.first;
  }

  inline bool operator()(const std::pair<Key, Value>& lhs, const Key& rhs) const
  {
    return lhs.first < rhs;
  }

  inline bool operator()(const Key& lhs, const std::pair<Key, Value>& rhs) const
  {
    return lhs < rhs.first;
  }

  inline bool operator()(const Key& lhs, const Key& rhs) const
  {
    return lhs < rhs;
  }
};

template<typename T>
class SortedKeyBoolPairVector {
public:
  SortedKeyBoolPairVector() = default;
  ~SortedKeyBoolPairVector() = default;

  bool set_if_not_already_set(const T& key, bool value)
  {
    auto iter = std::lower_bound(m_values.begin(), m_values.end(), key, m_compare);

    if (iter == m_values.end() || iter->first != key) {
      m_values.insert(iter, std::make_pair(key, value));
      return true;
    }

    return false;
  }

  bool set(const T& key, bool value)
  {
    auto iter = std::lower_bound(m_values.begin(), m_values.end(), key, m_compare);

    if (iter == m_values.end() || iter->first != key) {
      m_values.insert(iter, std::make_pair(key, value));
    }

    iter->second = value;
    return true;
  }

  bool get(const T& key) const
  {
    auto iter = std::lower_bound(m_values.begin(), m_values.end(), key, m_compare);

    if (iter == m_values.end() || iter->first != key) {
      return false;
    }

    return iter->second;
  }

  void clear()
  {
    m_values.clear();
  }

private:
  std::vector<std::pair<T, bool>> m_values;
  ComparePairByKey<T, bool> m_compare;
};

// Implementation of above using maps ... which is better to use?
template<typename T>
class SortedKeyBoolPairMap {
public:
  SortedKeyBoolPairMap() = default;
  ~SortedKeyBoolPairMap() = default;

  bool set_if_not_already_set(const T& key, bool value)
  {
    if(m_values.find(key) == m_values.end()) {
      m_values[key] = value;
      return true;
    }

    return false;
  }

  bool set(const T& key, bool value)
  {
    m_values[key] = value;
    return true;
  }

  bool get(const T& key)
  {
    return m_values[key];
  }

  void clear()
  {
    m_values.clear();
  }
private:
  std::map<T, bool> m_values;
};

////////////////////////////////////////////////////////////

template<typename T>
class SortedKeyIntVectorPairVector {
public:
  SortedKeyIntVectorPairVector() = default;
  ~SortedKeyIntVectorPairVector() = default;

  bool add(const T& key, int value)
  {
    auto iter = std::lower_bound(m_values.begin(), m_values.end(), key, m_compare);

    if (iter == m_values.end() || iter->first != key) {
      m_values.insert(iter, std::make_pair(key, std::vector<int>{value}));
      return true;
    }

    return stk::util::insert_keep_sorted_and_unique(value, iter->second);
  }

  const std::vector<int>& get(const T& key) const
  {
    auto iter = std::lower_bound(m_values.begin(), m_values.end(), key, m_compare);

    if (iter == m_values.end() || iter->first != key) {
      return m_emptyVector;
    }

    return iter->second;
  }

  void clear()
  {
    m_values.clear();
  }

private:
  std::vector<int> m_emptyVector;
  std::vector<std::pair<T, std::vector<int>>> m_values;
  ComparePairByKey<T, std::vector<int>> m_compare;
};

}}} // end namespaces stk mesh

#endif
