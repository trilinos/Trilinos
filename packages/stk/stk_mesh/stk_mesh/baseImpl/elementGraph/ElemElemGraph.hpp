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

#ifndef STK_ELEM_ELEM_GRAPH_HPP
#define STK_ELEM_ELEM_GRAPH_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/baseImpl/elementGraph/ParallelInfoForGraph.hpp>
#include <stk_mesh/base/SideIdPool.hpp>
#include <stk_mesh/baseImpl/DeletedElementInfo.hpp>
#include <stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp>

#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "SideConnector.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "BulkDataIdMapper.hpp"

namespace stk { class CommBuffer; }

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh {

struct moved_parallel_graph_info {

    moved_parallel_graph_info(int in_proc_to_tell,
                              stk::mesh::EntityId in_elem_id,
                              int in_elem_side,
                              stk::mesh::EntityId in_moved_elem_id,
                              int in_moved_elem_side,
                              int in_destination_proc)
     :
            elem_id(in_elem_id),
            moved_elem_id(in_moved_elem_id),
            elem_side(in_elem_side),
            moved_elem_side(in_moved_elem_side),
            proc_to_tell(in_proc_to_tell),
            destination_proc(in_destination_proc)
    {
    }

    stk::mesh::EntityId elem_id;
    stk::mesh::EntityId moved_elem_id;
    int elem_side;
    int moved_elem_side;
    int proc_to_tell;
    int destination_proc;
};

struct RemoteEdge;

typedef std::map<EntityVector, stk::mesh::impl::ParallelElementDataVector> SideNodeToReceivedElementDataMap;

class SharedSidesCommunication
{
public:
    SharedSidesCommunication(stk::mesh::BulkData& bulkData,
                             const impl::ElemSideProcVector & elementSidesToSend)
    : m_bulkData(bulkData),
      m_commSparse(bulkData.parallel()),
      m_elementSidesToSend(elementSidesToSend)
    {
    }

    void communicate_element_sides();
    void unpack_side_data_map(SideNodeToReceivedElementDataMap& elementSidesReceived);
    void unpack_side_data_vec(impl::ParallelElementDataVector& elementSidesReceived);

private:
    void pack_shared_side_nodes_of_elements(stk::CommSparse& comm) const;

    stk::mesh::BulkData& m_bulkData;
    stk::CommSparse m_commSparse;
    const impl::ElemSideProcVector & m_elementSidesToSend;
};


class ElemElemGraph
{
public:
    ElemElemGraph(stk::mesh::BulkData& bulkData);

    virtual ~ElemElemGraph();

    size_t get_num_connected_elems(stk::mesh::Entity localElement) const;

    bool is_connected_elem_locally_owned(stk::mesh::Entity localElement, size_t indexConnElement) const;

    impl::ElementViaSidePair get_connected_element_and_via_side(stk::mesh::Entity localElement, size_t indexConnElement) const;

    impl::IdViaSidePair get_connected_remote_id_and_via_side(stk::mesh::Entity localElement, size_t indexConnElement) const;

    int get_connected_elements_side(stk::mesh::Entity localElement, size_t indexConnElement) const;

    int get_owning_proc_id_of_remote_element(stk::mesh::Entity localElement, size_t indexConnElement) const;

    bool is_connected_to_other_element_via_side_ordinal(stk::mesh::Entity element, int sideOrdinal) const;

    impl::ParallelInfo& get_parallel_edge_info(stk::mesh::Entity element, int side1, stk::mesh::EntityId remote_id, int side2);
    const impl::ParallelInfo& get_const_parallel_edge_info(stk::mesh::Entity element, int side1, stk::mesh::EntityId remote_id, int side2) const;

    size_t num_edges() const;

    size_t num_parallel_edges() const;

    void add_elements(const stk::mesh::EntityVector &elements);

    void delete_elements(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete);

    size_t size() const {return m_graph.get_num_elements_in_graph() - m_deleted_element_local_id_pool.size();}

    impl::LocalId get_local_element_id(stk::mesh::Entity local_element, bool require_valid_id = true) const;

    void fill_from_mesh();

    stk::mesh::SideConnector& get_side_connector() { return m_sideConnector; }
    stk::mesh::SideNodeConnector& get_side_node_connector() { return m_sideNodeConnector; }
    stk::mesh::SideIdChooser get_side_id_chooser();

    const stk::mesh::BulkData& get_mesh() const;

    GraphEdgesForElement get_edges_for_element(impl::LocalId elem) const
    {
        return m_graph.get_edges_for_element(elem);
    }

    const std::vector<GraphEdge> & get_coincident_edges_for_element(impl::LocalId elem) const
    {
        return m_coincidentGraph.get_edges_for_element(elem);
    }

    std::vector<stk::mesh::impl::LocalId> get_coincident_elements(impl::LocalId elem) const
    {
        std::vector<stk::mesh::impl::LocalId> coincidents;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdges = get_coincident_edges_for_element(elem);
        for(const stk::mesh::GraphEdge& coinEdge : coincidentEdges)
            coincidents.push_back(coinEdge.elem2());
        stk::util::sort_and_unique(coincidents);
        return coincidents;
    }

    std::vector<stk::mesh::impl::LocalId> get_coincident_elements_on_side(impl::LocalId elem, int side) const
    {
        std::vector<stk::mesh::impl::LocalId> coincidents;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdges = get_coincident_edges_for_element(elem);
        for(const stk::mesh::GraphEdge& coinEdge : coincidentEdges) {
          if(coinEdge.side1() == side) {
            coincidents.push_back(coinEdge.elem2());
          }
        }
        stk::util::sort_and_unique(coincidents);
        return coincidents;
    }

    stk::mesh::EntityId convert_negative_local_id_to_global_id(impl::LocalId localId) const
    {
        return m_parallelInfoForGraphEdges.convert_negative_local_id_to_remote_global_id(localId);
    }

    const impl::ParallelInfo & get_parallel_info_for_graph_edge(const stk::mesh::GraphEdge& edge) const
    {
        return m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(edge);
    }

    const ParallelInfoForGraphEdges& get_parallel_info_for_graph_edges() const
    {
        return m_parallelInfoForGraphEdges;
    }

    unsigned mod_cycle_when_graph_modified() const { return m_modCycleWhenGraphModified; }

    const Graph& get_graph() const { return m_graph; }
    Graph& get_graph() { return m_graph; }

    stk::mesh::ParallelInfoForGraphEdges& get_parallel_graph() { return m_parallelInfoForGraphEdges; }
    const stk::mesh::ParallelInfoForGraphEdges& get_parallel_graph() const { return m_parallelInfoForGraphEdges; }

    void create_side_entities(const std::vector<int> &exposedSides,
                              impl::LocalId local_id,
                              const stk::mesh::PartVector& skin_parts,
                              std::vector<stk::mesh::sharing_info> &shared_modified);
    void write_graph(std::ostream& out, const std::string preamble = "") const;

    stk::mesh::Entity get_entity(stk::mesh::impl::LocalId localId) const
    {
        return m_idMapper.local_to_entity(localId);
    }

protected:
    void fill_graph();
    void fill_parallel_graph(impl::ElemSideProcVector & elementSidesToSend, impl::ParallelElementDataVector & elementSidesReceived);
    void create_parallel_graph_edge(const impl::ParallelElementData &elementData,
                                    const stk::mesh::impl::ParallelElementData &remoteElementData,
                                    impl::ElemSideProcVector & elementSidesToSend,
                                    std::vector<impl::SharedEdgeInfo> &newlySharedEdges,
                                    impl::ParallelGraphInfo& newParallelGraphInfo);

    void communicate_shared_sides(impl::ElemSideProcVector& elementSidesToSend,
                                  impl::ParallelElementDataVector& elementSidesReceived);

    stk::topology get_topology_of_connected_element(const GraphEdge &graphEdge);

    stk::topology get_topology_of_remote_element(const GraphEdge &graphEdge);

    void add_local_graph_edges_for_elem(const stk::mesh::MeshIndex &meshIndex, impl::LocalId local_elem_id, std::vector<stk::mesh::GraphEdge> &graphEdges,
                                        std::vector<stk::mesh::GraphEdge> &coincidentGraphEdges,
                                        stk::mesh::EntityVector& side_nodes,
                                        impl::SerialElementDataVector& connectedElementDataVector,
                                        bool only_consider_upper_symmetry = true) const;

    void fill_elements_attached_to_local_nodes(const stk::mesh::EntityVector& sideNodesOfReceivedElement,
                                               stk::mesh::EntityId elementId,
                                               stk::topology elementTopology,
                                               impl::SerialElementDataVector& connectedElementDataVector) const;
    void get_elements_attached_to_remote_nodes(const stk::mesh::EntityVector& sideNodesOfReceivedElement,
                                                                          stk::mesh::EntityId elementId,
                                                                          stk::topology elementTopology,
                                                                          impl::ParallelElementDataVector& connectedElementDataVector) const;

    impl::LocalId get_new_local_element_id_from_pool();
    int size_data_members();
    void clear_data_members();
    void pack_deleted_element_comm(stk::CommSparse &comm,
                                   const std::vector<impl::DeletedElementData> &local_elem_and_remote_connected_elem);

    void pack_shell_connectivity(stk::CommSparse & comm, const std::vector<impl::ShellConnectivityData> & shellConnectivityList,
                                 const std::vector<impl::ElementSidePair> &deletedShells);

    void pack_remote_edge_across_shell(stk::CommSparse &comm, stk::mesh::EntityVector &addedShells, int phase);

    void communicate_remote_edges_for_pre_existing_graph_nodes(const std::vector<impl::SharedEdgeInfo> &newlySharedEdges,
                                                          std::vector<impl::SharedEdgeInfo> &receivedSharedEdges);

    void connect_remote_element_to_existing_graph(const impl::SharedEdgeInfo &receivedSharedEdge);

    void collect_local_shell_connectivity_data(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete,
                                               std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                               std::vector<impl::ElementSidePair> &deletedShells);

    bool communicate_if_shell_connectivity(std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                            const std::vector<impl::ElementSidePair> &deletedShells);

    void delete_local_connections_and_collect_remote(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete,
                                                     std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem);

    void communicate_remote_connections_to_delete(const std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem,
                                                  std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges);

    void clear_deleted_element_connections(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete);

    void delete_remote_connections(const std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges);

    void delete_remote_connection(stk::mesh::Entity connected_elem,
                                  stk::mesh::EntityId deleted_elem_global_id,
                                  std::vector<GraphEdge>& pllEdgesToDelete);

    void reconnect_volume_elements_across_deleted_shells(std::vector<impl::ShellConnectivityData> & shellConnectivityList);

    void break_remote_shell_connectivity_and_pack(stk::CommSparse& comm, impl::LocalId leftId, impl::LocalId rightId, int phase);

    void pack_both_remote_shell_connectivity(stk::CommSparse &comm, impl::LocalId shellId, impl::LocalId leftId, impl::LocalId rightId);

    void unpack_remote_edge_across_shell(stk::CommSparse &comm);

    bool is_valid_graph_element(stk::mesh::Entity local_element) const;

    bool is_valid_graph_edge(const GraphEdge &graphEdge) const;

    stk::mesh::BulkData &m_bulk_data;
    Graph m_graph;
    std::vector<GraphEdge> m_edgesToAdd;
    std::vector<GraphEdge> m_edgesToDelete;
    std::vector<GraphEdge> m_newGraphEdges;
    std::vector<GraphEdge> m_coincidentGraphEdges;
    ParallelInfoForGraphEdges m_parallelInfoForGraphEdges;
    unsigned m_modCycleWhenGraphModified;
    std::vector<impl::LocalId> m_deleted_element_local_id_pool;
    std::vector<stk::topology> m_element_topologies;
    bool m_any_shell_elements_exist;
    std::vector<int> m_deleted_elem_pool;
    mutable stk::mesh::EntityVector m_scratchEntityVector;
    impl::SparseGraph m_coincidentGraph;
    impl::ElementLocalIdMapper m_idMapper;
    SideConnector m_sideConnector;
    SideNodeConnector m_sideNodeConnector;
private:
    void add_side_for_remote_edge(const GraphEdge & graphEdge,
                                                 int elemSide,
                                                 stk::mesh::Entity element,
                                                 const stk::mesh::PartVector& skin_parts,
                                                 std::vector<stk::mesh::sharing_info> &shared_modified);

    template <typename GraphType>
    void add_new_local_edges_to_graph(GraphType &graph, const std::vector<stk::mesh::GraphEdge> &newGraphEdges);

    void add_vertex(impl::LocalId new_elem_id, stk::mesh::Entity elem_to_add);
    std::pair<stk::mesh::EntityVector,stk::mesh::EntityVector> filter_add_elements_arguments(const stk::mesh::EntityVector& allUnfilteredElementsNotAlreadyInGraph) const;
    stk::mesh::impl::DeletedElementInfoVector filter_delete_elements_argument(const stk::mesh::impl::DeletedElementInfoVector& elements_to_delete_argument) const;
    void add_elements_locally(const stk::mesh::EntityVector& allElementsNotAlreadyInGraph);
    stk::mesh::Entity add_side_to_mesh(const stk::mesh::impl::ElementSidePair& side_pair, const stk::mesh::PartVector& skin_parts);

    void write_graph_edge(std::ostringstream& os, const stk::mesh::GraphEdge& graphEdge) const;
    unsigned get_max_num_sides_per_element() const;
    bool did_already_delete_a_shell_between_these_elements(std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                                           const impl::ShellConnectivityData& shellConnectivityData);
    bool are_connectivities_for_same_graph_edge(const impl::ShellConnectivityData& shellConn,
                                                      const impl::ShellConnectivityData& shellConnectivityData);
    bool is_connected_to_shell_on_side(stk::mesh::impl::LocalId localElemLocalId, int side);
    void insert_edge_between_elements(impl::LocalId local_elem_id,
                                      int side_index,
                                      const impl::SerialElementData& elemData,
                                      std::vector<stk::mesh::GraphEdge>& coincidentGraphEdges) const;

    std::string print_edge(const GraphEdge& graphEdge);
};

namespace impl
{

template<typename SideData>
void add_second_side_if_shell(const stk::mesh::Entity* connectedElemNodes, SideData& elemData, std::vector<SideData> & element_data)
{
    if (elemData.get_element_topology().is_shell())
    {
        // Also add the back-side face if this is a shell
        unsigned index = (elemData.get_element_side_index() == 0) ? 1 : 0;
        elemData.set_element_side_index(index);
        elemData.get_element_topology().side_nodes(connectedElemNodes, elemData.get_element_side_index(), elemData.side_nodes_begin());
        element_data.push_back(elemData);
    }
}

inline
stk::mesh::OrdinalAndPermutation flip_shell_to_get_opposing_normal(const stk::mesh::OrdinalAndPermutation &connectedOrdAndPerm,
                                                                        stk::topology topology)
{
    STK_ThrowRequireWithSierraHelpMsg(static_cast<unsigned>(connectedOrdAndPerm.second) < topology.num_positive_permutations());
    unsigned sideOrdinal = connectedOrdAndPerm.first == 0u ? 1u : 0u;
    unsigned perm = connectedOrdAndPerm.second + topology.num_positive_permutations();
    return stk::mesh::OrdinalAndPermutation(static_cast<stk::mesh::ConnectivityOrdinal>(sideOrdinal),
                                            static_cast<stk::mesh::Permutation>(perm));
}

impl::ElemSideProcVector gather_element_side_ids_to_send(const stk::mesh::BulkData& bulkData);

stk::mesh::EntityVector gather_solid_elements_connected_to_shell(const stk::mesh::BulkData& bulkData, stk::mesh::Entity shellElement);

} // end impl

}} // end stk mesh namespaces

#endif
