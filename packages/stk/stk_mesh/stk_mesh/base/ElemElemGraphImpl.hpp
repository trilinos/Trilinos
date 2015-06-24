#ifndef STK_ELEM_ELEM_GRAPH_IMPL_HPP
#define STK_ELEM_ELEM_GRAPH_IMPL_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { struct sharing_info; } }
namespace stk { namespace mesh { class ElemElemGraph; } }
namespace stk { class CommSparse; }

namespace stk { namespace mesh {
namespace impl
{

typedef int64_t LocalId;

struct graphEdgeProc
{
    stk::mesh::EntityId m_localElementId;
    stk::mesh::EntityId m_remoteElementId;
    int m_proc_id;
    graphEdgeProc(const stk::mesh::EntityId& localElementId, const stk::mesh::EntityId &remoteElementId, int proc_id) :
        m_localElementId(localElementId), m_remoteElementId(remoteElementId), m_proc_id(proc_id) {}
};

struct parallel_info
{
    int m_other_proc;
    int m_other_side_ord;
    int m_permutation;
    bool m_in_part;
    stk::mesh::EntityId m_chosen_side_id;

    parallel_info(int proc, int side_ord, int perm, stk::mesh::EntityId chosen_face_id) :
        m_other_proc(proc), m_other_side_ord(side_ord), m_permutation(perm), m_in_part(true),
        m_chosen_side_id(chosen_face_id) {}
};

struct ConnectedElementData
{
    int m_procId;
    LocalId m_elementId;
    stk::topology m_elementTopology;
    unsigned m_sideIndex;
    stk::mesh::EntityId m_suggestedFaceId;
    stk::mesh::EntityVector m_sideNodes;
};

typedef std::pair<LocalId,int> ElementSidePair;
typedef std::map<std::pair<LocalId,stk::mesh::EntityId>, parallel_info > ParallelGraphInfo;
typedef std::vector<std::vector<LocalId> > ElementGraph;
typedef std::vector<std::vector<int> > SidesForElementGraph;
typedef std::vector<ConnectedElementData> ConnectedElementDataVector;


NAMED_PAIR( EntitySidePair , stk::mesh::Entity , entity , unsigned , side_id )
NAMED_PAIR( ProcFaceIdPair , int , proc , stk::mesh::EntityId , side_id )

typedef std::multimap<EntitySidePair, ProcFaceIdPair>  ElemSideToProcAndFaceId;

void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<stk::topology>& element_topologies);
void fill_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<LocalId>& entity_to_local_id, std::vector<stk::topology>& element_topologies);

ElemSideToProcAndFaceId get_element_side_ids_to_communicate(const stk::mesh::BulkData& bulkData);

void pack_shared_side_nodes_of_elements(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, ElemSideToProcAndFaceId& elements_to_communicate,
        const std::vector<stk::mesh::EntityId>& suggested_face_ids);

void add_possibly_connected_elements_to_graph_using_side_nodes(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, ParallelGraphInfo& parallel_graph_info,
        const ElemSideToProcAndFaceId& elemSideComm,
        std::vector<ConnectedElementData> & communicatedElementDataVector);

void fill_parallel_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, ParallelGraphInfo& parallel_graph_info,
        ElemSideToProcAndFaceId& elem_side_comm, const std::vector<stk::mesh::EntityId>& suggested_face_ids);

std::vector<graphEdgeProc> get_elements_to_communicate(stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &killedElements,
        const ElemElemGraph& elem_graph);

void communicate_killed_entities(stk::mesh::BulkData& bulkData, const std::vector<graphEdgeProc>& elements_to_comm,
        std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges);

void pack_elements_to_comm(stk::CommSparse &comm, const std::vector<graphEdgeProc>& elements_to_comm);

bool create_or_delete_shared_side(stk::mesh::BulkData& bulkData, const parallel_info& parallel_edge_info, const ElemElemGraph& elementGraph,
        stk::mesh::Entity local_element, stk::mesh::EntityId remote_id, bool create_face, const stk::mesh::PartVector& face_parts,
        std::vector<stk::mesh::sharing_info> &shared_modified, stk::mesh::EntityVector &deletedEntities,
        size_t &id_counter, stk::mesh::EntityId suggested_local_face_id, stk::mesh::Part& faces_created_during_death);

stk::mesh::Entity get_side_for_element(const stk::mesh::BulkData& bulkData, stk::mesh::Entity this_elem_entity, int side_id);

int get_element_side_multiplier();

bool is_id_already_in_use_locally(stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank, stk::mesh::EntityId id);

bool does_side_exist_with_different_permutation(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::ConnectivityOrdinal side_ordinal, stk::mesh::Permutation perm);

bool does_element_side_exist(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal);

stk::mesh::Entity connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
        stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts);

stk::mesh::EntityId get_side_global_id(const stk::mesh::BulkData &bulkData, const ElemElemGraph& elementGraph, stk::mesh::Entity element1, stk::mesh::Entity element2,
        int element1_side_id);

void filter_for_candidate_elements_to_connect(const stk::mesh::BulkData & mesh,
                                          stk::mesh::Entity localElement,
                                          const unsigned sideOrdinal,
                                          ConnectedElementDataVector & connectedElementData);

void break_volume_element_connections_across_shells(const std::set<EntityId> & localElementsConnectedToRemoteShell,
                                       ElementGraph & elem_graph,
                                       SidesForElementGraph & via_sides);
}
}} // end namespaces stk mesh

#endif
