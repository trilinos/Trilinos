#ifndef STK_ELEM_ELEM_GRAPH_HPP
#define STK_ELEM_ELEM_GRAPH_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SideIdPool.hpp>
#include <stk_mesh/baseImpl/DeletedElementInfo.hpp>
#include <stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp>

#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "GraphEdgeData.hpp"
#include "SideConnector.hpp"
#include "../MeshImplUtils.hpp"
#include "../../../../stk_util/stk_util/util/SortAndUnique.hpp"

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

void change_entity_owner(stk::mesh::BulkData &bulkData, stk::mesh::ElemElemGraph &elem_graph,
                         std::vector< std::pair< stk::mesh::Entity, int > > &elem_proc_pairs_to_move,
                         stk::mesh::Part *active_part=NULL);

struct RemoteEdge;

class ElemElemGraph
{
public:
    ElemElemGraph(stk::mesh::BulkData& bulkData, const stk::mesh::Selector &selector, const stk::mesh::Selector *air = nullptr);

    virtual ~ElemElemGraph();

    size_t get_num_connected_elems(stk::mesh::Entity localElement) const;

    unsigned get_num_connected_elems(stk::mesh::Entity localElement, int side_id) const;

    bool is_connected_elem_locally_owned(stk::mesh::Entity localElement, size_t indexConnElement) const;

    impl::ElementViaSidePair get_connected_element_and_via_side(stk::mesh::Entity localElement, size_t indexConnElement) const;

    impl::IdViaSidePair get_connected_remote_id_and_via_side(stk::mesh::Entity localElement, size_t indexConnElement) const;

    int get_connected_elements_side(stk::mesh::Entity localElement, size_t indexConnElement) const;

    int get_owning_proc_id_of_remote_element(stk::mesh::Entity localElement, size_t indexConnElement) const;

    bool is_connected_to_other_element_via_side_ordinal(stk::mesh::Entity element, int sideOrdinal) const;

    impl::ParallelInfo& get_parallel_edge_info(stk::mesh::Entity element, int side1, stk::mesh::EntityId remote_id, int side2);
    const impl::ParallelInfo& get_const_parallel_edge_info(stk::mesh::Entity element, int side1, stk::mesh::EntityId remote_id, int side2) const;

    size_t num_edges() const;

    size_t num_parallel_edges() const { return m_num_parallel_edges; }

    void add_elements(const stk::mesh::EntityVector &elements);

    void delete_elements(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete);

    size_t size() {return m_graph.get_num_elements_in_graph() - m_deleted_element_local_id_pool.size();}

    impl::LocalId get_local_element_id(stk::mesh::Entity local_element, bool require_valid_id = true) const;

    void fill_from_mesh();
    void create_parallel_graph_info_needed_once_entities_are_moved(const stk::mesh::EntityProcVec &elemProcPairsToMove,
                                         impl::ParallelGraphInfo &new_parallel_graph_entries);
    stk::mesh::EntityId get_available_side_id();

    stk::mesh::SideConnector get_side_connector();

    const stk::mesh::BulkData& get_mesh() const;

    const GraphEdgesForElement& get_edges_for_element(impl::LocalId elem) const
    {
        return m_graph.get_edges_for_element(elem);
    }

    const std::vector<GraphEdge> & get_coincident_edges_for_element(impl::LocalId elem) const
    {
        return m_coincidentGraph.get_edges_for_element(elem);
    }

    stk::mesh::Entity get_entity_from_local_id(impl::LocalId elem) const
    {
        return m_local_id_to_element_entity[elem];
    }

    const impl::ParallelInfo & get_parallel_info_for_graph_edge(const stk::mesh::GraphEdge& edge) const
    {
        return m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(edge);
    }

    const ParallelInfoForGraphEdges& get_parallel_info_for_graph_edges() const
    {
        return m_parallelInfoForGraphEdges;
    }

    const Graph& get_graph() const
    {
        return m_graph;
    }

    const stk::mesh::ParallelInfoForGraphEdges& get_parallel_graph() const { return m_parallelInfoForGraphEdges; }

    void create_side_entities(const std::vector<int> &exposedSides,
                              impl::LocalId local_id,
                              const stk::mesh::PartVector& skin_parts,
                              std::vector<stk::mesh::sharing_info> &shared_modified);
    void write_graph(std::ostream& out) const;

    stk::mesh::Entity get_entity(stk::mesh::impl::LocalId local_id) const { return m_local_id_to_element_entity[local_id]; }

    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > get_split_coincident_elements();
protected:
    void fill_graph();
    void update_number_of_parallel_edges();
    void fill_parallel_graph(impl::ElemSideToProcAndFaceId& elem_side_comm);

    void add_possibly_connected_elements_to_graph_using_side_nodes( const stk::mesh::impl::ElemSideToProcAndFaceId& elemSideComm,
                                                                    stk::mesh::impl::ParallelElementDataVector & communicatedElementDataVector,
                                                                    std::vector<impl::SharedEdgeInfo> &newlySharedEdges);

    stk::mesh::EntityId pick_id_for_side_if_created(const impl::ParallelElementData & elemDataFromOtherProc,
            const stk::mesh::impl::ElemSideToProcAndFaceId& elemSideDataSent, stk::mesh::Entity localElem, unsigned side_index);

    void add_parallel_edge_and_info(const stk::mesh::impl::ElemSideToProcAndFaceId& elemSideDataSent,
            const impl::ParallelElementDataVector &filteredCommunicatedElementData, const impl::ParallelElementData &elementData,
            std::vector<impl::SharedEdgeInfo> &newlySharedEdges);

    stk::topology get_topology_of_connected_element(const GraphEdge &graphEdge);

    stk::topology get_topology_of_remote_element(const GraphEdge &graphEdge);

    void add_local_graph_edges_for_elem(const stk::mesh::MeshIndex &meshIndex, impl::LocalId local_elem_id, std::vector<stk::mesh::GraphEdge> &graphEdges) const;

    impl::SerialElementDataVector get_only_valid_element_connections(stk::mesh::Entity element, unsigned side_index, const stk::mesh::EntityVector& side_nodes) const;

    impl::LocalId create_new_local_id(stk::mesh::Entity new_elem);

    impl::LocalId get_new_local_element_id_from_pool();
    int size_data_members();
    void clear_data_members();
    void resize_entity_to_local_id_if_needed(size_t max_index);
    size_t find_max_local_offset_in_neighborhood(stk::mesh::Entity element);
    void pack_deleted_element_comm(stk::CommSparse &comm,
                                   const std::vector<impl::DeletedElementData> &local_elem_and_remote_connected_elem);

    void pack_shell_connectivity(stk::CommSparse & comm, const std::vector<impl::ShellConnectivityData> & shellConnectivityList,
                                 const std::vector<impl::ElementSidePair> &deletedShells);

    void pack_remote_edge_across_shell(stk::CommSparse &comm, stk::mesh::EntityVector &addedShells, int phase);

    void filter_for_elements_in_graph(stk::mesh::EntityVector &localElements);

    void communicate_remote_edges_for_pre_existing_graph_nodes(const std::vector<impl::SharedEdgeInfo> &newlySharedEdges,
                                                          std::vector<impl::SharedEdgeInfo> &receivedSharedEdges);

    void connect_remote_element_to_existing_graph(const impl::SharedEdgeInfo &receivedSharedEdge);

    void collect_local_shell_connectivity_data(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete,
                                               std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                               std::vector<impl::ElementSidePair> &deletedShells);

    void communicate_shell_connectivity(std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                        const std::vector<impl::ElementSidePair> &deletedShells);

    void delete_local_connections_and_collect_remote(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete,
                                                     std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem);

    void communicate_remote_connections_to_delete(const std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem,
                                                  std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges);

    void clear_deleted_element_connections(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete);

    void delete_remote_connections(const std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges);

    void reconnect_volume_elements_across_deleted_shells(std::vector<impl::ShellConnectivityData> & shellConnectivityList);

    void break_remote_shell_connectivity_and_pack(stk::CommSparse& comm, impl::LocalId leftId, impl::LocalId rightId, int phase);

    void pack_both_remote_shell_connectivity(stk::CommSparse &comm, impl::LocalId shellId, impl::LocalId leftId, impl::LocalId rightId);

    void unpack_remote_edge_across_shell(stk::CommSparse &comm);

    bool is_valid_graph_element(stk::mesh::Entity local_element) const;

    stk::mesh::BulkData &m_bulk_data;
    const stk::mesh::Selector m_skinned_selector;
    const stk::mesh::Selector* m_air_selector;
    Graph m_graph;
    ParallelInfoForGraphEdges m_parallelInfoForGraphEdges;
    stk::mesh::EntityVector m_local_id_to_element_entity;
    std::vector<impl::LocalId> m_entity_to_local_id;
    std::vector<impl::LocalId> m_deleted_element_local_id_pool;
    std::vector<bool> m_local_id_in_pool;
    std::vector<stk::topology> m_element_topologies;
    std::vector<int> m_deleted_elem_pool;
    size_t m_num_parallel_edges;
    stk::mesh::SideIdPool m_sideIdPool;
    impl::SparseGraph m_coincidentGraph;
    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > m_splitCoincidents;
private:
    void resize_entity_to_local_id_vector_for_new_elements(const stk::mesh::EntityVector& allElementsNotAlreadyInGraph);
    stk::mesh::EntityId add_side_for_remote_edge(const GraphEdge & graphEdge,
                                                 int elemSide,
                                                 stk::mesh::Entity element,
                                                 const stk::mesh::PartVector& skin_parts,
                                                 std::vector<stk::mesh::sharing_info> &shared_modified);
    void add_both_edges_between_local_elements(const GraphEdge& graphEdge);
    void add_local_edges(stk::mesh::Entity elem_to_add, impl::LocalId new_elem_id);
    void add_vertex(impl::LocalId new_elem_id, stk::mesh::Entity elem_to_add);
    stk::mesh::EntityVector filter_add_elements_arguments(const stk::mesh::EntityVector& allUnfilteredElementsNotAlreadyInGraph) const;
    stk::mesh::impl::DeletedElementInfoVector filter_delete_elements_argument(const stk::mesh::impl::DeletedElementInfoVector& elements_to_delete_argument) const;
    impl::ElemSideToProcAndFaceId get_element_side_ids_to_communicate() const;
    void add_elements_locally(const stk::mesh::EntityVector& allElementsNotAlreadyInGraph);
    stk::mesh::Entity add_side_to_mesh(const stk::mesh::impl::ElementSidePair& side_pair, const stk::mesh::PartVector& skin_parts);

    void create_remote_sides(stk::mesh::BulkData& bulk_data,
                             const std::vector<RemoteEdge>& remote_edges,
                             stk::mesh::EntityVector& skinned_elements,
                             const stk::mesh::PartVector& skin_parts,
                             const std::vector<unsigned>& side_counts,
                             std::vector<stk::mesh::sharing_info>& shared_modified);
    void create_remote_sides1(stk::mesh::BulkData& bulk_data,
                              const std::vector<RemoteEdge>& remote_edges,
                              stk::mesh::EntityVector& skinned_elements,
                              const stk::mesh::PartVector& skin_parts,
                              const std::vector<unsigned>& side_counts,
                              std::vector<stk::mesh::sharing_info>& shared_modified);

    void extract_split_coincident_connections(stk::mesh::Entity localElem, const impl::ParallelElementData &localElementData, const impl::ParallelElementDataVector &elementsConnectedToThisElementSide);
    void add_split_coincident_connection(stk::mesh::Entity localElem, const impl::ParallelElementData & connectedElemParallelData);

    void update_all_local_neighbors(const stk::mesh::Entity elemToSend,
                                    const int destination_proc,
                                    impl::ParallelGraphInfo &newParallelGraphEntries);

    impl::ParallelInfo create_parallel_info(stk::mesh::Entity connected_element,
                                                            const stk::mesh::Entity elemToSend,
                                                            int side_id,
                                                            const int destination_proc);

    stk::mesh::Permutation get_permutation_given_neighbors_node_ordering(stk::mesh::Entity connected_element,
                                                                                        const stk::mesh::Entity elemToSend,
                                                                                        int side_id);

    size_t pack_shared_side_nodes_of_elements(stk::CommSparse &comm, impl::ElemSideToProcAndFaceId& elements_to_communicate);
    void extract_coincident_edges_and_fix_chosen_side_ids();
    std::vector<impl::LocalId> get_local_ids_for_element_entities(const stk::mesh::EntityVector &elems);

    void write_graph_edge(std::ostringstream& os, const stk::mesh::GraphEdge& graphEdge) const;
    unsigned get_max_num_sides_per_element() const;
    bool did_already_delete_a_shell_between_these_elements(std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                                           const impl::ShellConnectivityData& shellConnectivityData);
    bool are_connectivities_for_same_graph_edge(const impl::ShellConnectivityData& shellConn,
                                                      const impl::ShellConnectivityData& shellConnectivityData);
    bool is_connected_to_shell_on_side(stk::mesh::impl::LocalId localElemLocalId, int side);
};

bool process_killed_elements(stk::mesh::BulkData& bulkData, ElemElemGraph& elementGraph, const stk::mesh::EntityVector& killedElements, stk::mesh::Part& active,
        const stk::mesh::PartVector& side_parts, const stk::mesh::PartVector* boundary_mesh_parts = nullptr);

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

template<typename SideData>
void add_local_elements_to_connected_list(const stk::mesh::BulkData& bulkData,
                                          const stk::mesh::ElemElemGraph& graph,
                                          const stk::mesh::EntityVector & local_elements_attached_to_side_nodes,
                                          const stk::mesh::EntityVector & side_nodes_received,
                                          std::vector<SideData> & element_data)
{
    for (const stk::mesh::Entity element : local_elements_attached_to_side_nodes)
    {
        stk::mesh::OrdinalAndPermutation connectedOrdAndPerm = stk::mesh::get_ordinal_and_permutation(bulkData, element, bulkData.mesh_meta_data().side_rank(), side_nodes_received);

        if (INVALID_CONNECTIVITY_ORDINAL != connectedOrdAndPerm.first)
        {
            const stk::mesh::Bucket & connectedBucket = bulkData.bucket(element);
            const stk::mesh::Entity* connectedElemNodes = bulkData.begin_nodes(element);

            ThrowAssertMsg(connectedBucket.topology().side_topology(connectedOrdAndPerm.first).num_nodes() == side_nodes_received.size(),
                          "Error, number of nodes on sides of adjacent elements do not agree:  " <<
                           side_nodes_received.size() << " != " << connectedBucket.topology().side_topology(connectedOrdAndPerm.first).num_nodes());

            impl::LocalId local_id = graph.get_local_element_id(element, false);
            if (local_id != impl::INVALID_LOCAL_ID)
            {
                SideData elemData;
                elemData.set_element_local_id(local_id);
                elemData.set_element_identifier(bulkData.identifier(element));
                elemData.set_element_topology(connectedBucket.topology());
                elemData.set_element_side_index(connectedOrdAndPerm.first);
                elemData.resize_side_nodes(side_nodes_received.size());
                elemData.get_element_topology().side_nodes(connectedElemNodes, elemData.get_element_side_index(), elemData.side_nodes_begin());
                element_data.push_back(elemData);
                add_second_side_if_shell(connectedElemNodes, elemData, element_data);
            }
        }
    }
}

template<typename SideData>
std::vector<SideData> get_elements_connected_via_sidenodes(const stk::mesh::BulkData& bulk_data, const ElemElemGraph& graph, const stk::mesh::EntityVector &sideNodesOfReceivedElement)
{
    unsigned num_side_nodes_received_element = sideNodesOfReceivedElement.size();
    stk::mesh::EntityVector localElementsConnectedToReceivedSideNodes;
    impl::find_locally_owned_elements_these_nodes_have_in_common(bulk_data, num_side_nodes_received_element, sideNodesOfReceivedElement.data(), localElementsConnectedToReceivedSideNodes);
    stk::util::sort_and_unique(localElementsConnectedToReceivedSideNodes);

    std::vector<SideData> connectedElementDataVector;
    impl::add_local_elements_to_connected_list(bulk_data, graph, localElementsConnectedToReceivedSideNodes, sideNodesOfReceivedElement, connectedElementDataVector);

    return connectedElementDataVector;
}

} // end impl

}} // end stk mesh namespaces

#endif
