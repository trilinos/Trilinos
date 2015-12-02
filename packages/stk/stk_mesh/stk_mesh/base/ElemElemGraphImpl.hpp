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
namespace stk { namespace mesh { class Graph; } }
namespace stk { namespace mesh { class ParallelInfoForGraphEdges; } }
namespace stk { class CommSparse; }

namespace stk { namespace mesh {
namespace impl
{

typedef int64_t LocalId;

struct graphEdgeProc
{
    stk::mesh::EntityId m_localElementId;
    stk::mesh::EntityId m_remoteElementId;
    int m_localSide;
    int m_remoteSide;
    int m_proc_id;
    graphEdgeProc(const stk::mesh::EntityId& localElementId, int localSide, const stk::mesh::EntityId &remoteElementId, int remoteSide, int proc_id) :
        m_localElementId(localElementId), m_remoteElementId(remoteElementId), m_localSide(localSide), m_remoteSide(remoteSide), m_proc_id(proc_id) {}
};

struct parallel_info
{
    int m_other_proc;
    int m_permutation;
    stk::topology m_remote_element_toplogy;
    bool m_in_body_to_be_skinned;
    bool m_is_air;

    stk::mesh::EntityId m_chosen_side_id;

    parallel_info(int proc, int perm, stk::mesh::EntityId chosen_face_id, stk::topology other_elem_topology, bool inPart, bool isInAir=false) :
        m_other_proc(proc), m_permutation(perm), m_remote_element_toplogy(other_elem_topology), m_in_body_to_be_skinned(inPart), m_is_air(isInAir),
        m_chosen_side_id(chosen_face_id) {}
};

struct ConnectedElementData
{
    int m_procId;
    LocalId m_elementLocalId;
    stk::mesh::EntityId m_elementIdentifier;
    stk::topology m_elementTopology;
    unsigned m_sideIndex;
    stk::mesh::EntityId m_suggestedFaceId;
    stk::mesh::EntityVector m_sideNodes;
    bool m_isInPart;
    bool m_isAir;

    ConnectedElementData()
    : m_procId(-1),
      m_elementLocalId(std::numeric_limits<impl::LocalId>::max()),
      m_elementIdentifier(stk::mesh::InvalidEntityId),
      m_elementTopology(stk::topology::INVALID_TOPOLOGY),
      m_sideIndex(std::numeric_limits<unsigned>::max()),
      m_suggestedFaceId(stk::mesh::InvalidEntityId),
      m_isInPart(true), m_isAir(false)
    {}
};

struct SharedEdgeInfo
{
    stk::mesh::EntityId m_locaElementlId;
    stk::mesh::EntityId m_remoteElementId;
    int m_procId;
    unsigned m_sideIndex;
    stk::mesh::EntityId m_chosenSideId;
    stk::mesh::EntityVector m_sharedNodes;
    bool m_isInPart;
    bool m_isInAir;
    stk::topology m_remoteElementTopology;
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

struct GraphEdge
{
    GraphEdge(impl::LocalId e1, int s1, impl::LocalId e2, int s2) :
        elem1(e1), elem2(e2), side1(s1), side2(s2)
    {}
    GraphEdge() :
        elem1(std::numeric_limits<impl::LocalId>::max()), elem2(std::numeric_limits<impl::LocalId>::max()), side1(-1), side2(-1)
    {}
    impl::LocalId elem1;
    impl::LocalId elem2;
    int side1;
    int side2;
};

typedef GraphEdge CoincidentElementConnection;

struct GraphEdgeLessByElem2 {
    bool operator()(const GraphEdge& a, const GraphEdge& b) const
    {
        if (a.elem2 != b.elem2)
        {
            return a.elem2 < b.elem2;
        }
        else if (a.side2 != b.side2)
        {
            return a.side2 < b.side2;
        }
        else if (a.elem1 != b.elem1)
        {
            return a.elem1 < b.elem1;
        }
        else
        {
            return a.side1 < b.side1;
        }
    }
};

inline
bool operator==(const GraphEdge& a, const GraphEdge& b)
{
    return  a.elem1 == b.elem1 &&
            a.side1 == b.side1 &&
            a.elem2 == b.elem2 &&
            a.side2 == b.side2;
}

inline
std::ostream& operator<<(std::ostream& out, const GraphEdge& graphEdge)
{
    out << "(" << graphEdge.elem1 << "," << graphEdge.side1 << " -> " << graphEdge.elem2 << "," << graphEdge.side2 << ")";
    return out;
}

namespace impl {

typedef std::pair<LocalId,int> ElementSidePair;
typedef std::map<GraphEdge, parallel_info, GraphEdgeLessByElem2> ParallelGraphInfo;
typedef std::vector<std::vector<LocalId> > ElementGraph;
typedef std::vector<std::vector<int> > SidesForElementGraph;
typedef std::vector<ConnectedElementData> ConnectedElementDataVector;

NAMED_PAIR( EntitySidePair , stk::mesh::Entity , entity , unsigned , side_id )
NAMED_PAIR( ProcFaceIdPair , int , proc , stk::mesh::EntityId , side_id )
NAMED_PAIR( ProcVecFaceIdPair , std::vector<int> , proc_vec , stk::mesh::EntityId , side_id )

typedef std::multimap<EntitySidePair, ProcFaceIdPair>  ElemSideToProcAndFaceId;

void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<stk::topology>& element_topologies);
void fill_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<LocalId>& entity_to_local_id, std::vector<stk::topology>& element_topologies);

ElemSideToProcAndFaceId get_element_side_ids_to_communicate(const stk::mesh::BulkData& bulkData);
ElemSideToProcAndFaceId get_element_side_ids_to_communicate(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &element_list);

ElemSideToProcAndFaceId build_element_side_ids_to_proc_map(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &elements_to_communicate);

std::vector<graphEdgeProc> get_elements_to_communicate(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &killedElements,
        const ElemElemGraph& elem_graph);

std::vector<graphEdgeProc> communicate_killed_entities(stk::ParallelMachine communicator, const std::vector<graphEdgeProc>& elements_to_comm);

void pack_elements_to_comm(stk::CommSparse &comm, const std::vector<graphEdgeProc>& elements_to_comm);

void add_side_into_exposed_boundary(stk::mesh::BulkData& bulkData, const parallel_info& parallel_edge_info,
        stk::mesh::Entity local_element, int side_id, stk::mesh::EntityId remote_id, const stk::mesh::PartVector& parts_for_creating_side,
        std::vector<stk::mesh::sharing_info> &shared_modified, const stk::mesh::PartVector *boundary_mesh_parts = nullptr);

void remove_side_from_death_boundary(stk::mesh::BulkData& bulkData, stk::mesh::Entity local_element,
        stk::mesh::Part &activePart, stk::mesh::EntityVector &deletedEntities, int side_id);

stk::mesh::PartVector get_stk_parts_for_moving_parts_into_death_boundary(const stk::mesh::PartVector *bc_mesh_parts);

stk::mesh::Entity get_side_for_element(const stk::mesh::BulkData& bulkData, stk::mesh::Entity this_elem_entity, int side_id);

int get_element_side_multiplier();

bool is_id_already_in_use_locally(stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank, stk::mesh::EntityId id);

bool does_side_exist_with_different_permutation(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::ConnectivityOrdinal side_ordinal, stk::mesh::Permutation perm);

bool does_element_side_exist(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal);

stk::mesh::Entity connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
        stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts);

void filter_out_invalid_solid_shell_connections(const stk::mesh::BulkData & mesh,
                                          stk::mesh::Entity localElement,
                                          const unsigned sideOrdinal,
                                          ConnectedElementDataVector & connectedElementData);

void break_volume_element_connections_across_shells(const std::set<EntityId> & localElementsConnectedToRemoteShell,
                                       ElementGraph & elem_graph,
                                       SidesForElementGraph & via_sides);

void pack_newly_shared_remote_edges(stk::CommSparse &comm, const stk::mesh::BulkData &m_bulk_data, const std::vector<SharedEdgeInfo> &newlySharedEdges);

bool does_element_have_side(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element);

void add_element_side_pairs_for_unused_sides(LocalId elementId, size_t numElemSides, const stk::mesh::Graph& graph,
        std::vector<ElementSidePair>& element_side_pairs);

void create_sides_created_during_death_part(stk::mesh::MetaData &metaData);
stk::mesh::Part* get_sides_created_during_death_part(const stk::mesh::MetaData &metaData);
void add_parts_from_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::PartVector& side_parts);
stk::mesh::PartVector get_parts_for_creating_side(stk::mesh::BulkData& bulkData, const stk::mesh::PartVector& parts_for_creating_side, stk::mesh::Entity element, int side_ord);
bool side_created_during_death(stk::mesh::BulkData& bulkData, stk::mesh::Entity side);

struct CoincidentElementDescription
{
    int numSides;
    stk::mesh::impl::LocalId elem1;
    stk::mesh::impl::LocalId elem2;
};
typedef std::map<stk::mesh::impl::LocalId, std::vector<stk::mesh::GraphEdge>> SparseGraph;

bool are_elements_coincident(const stk::mesh::Graph &graph, const CoincidentElementDescription &elemDesc);
SparseGraph extract_coincident_sides(stk::mesh::Graph &graph, const std::vector<stk::topology> &topologies);
bool is_local_element(stk::mesh::impl::LocalId elemId);
std::map<int, std::vector<LocalId>> get_extracted_coincident_local_ids(const stk::mesh::Graph &graph,
                                                                       const stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                                                       const stk::mesh::impl::SparseGraph &extractedCoincidentElements);
void pack_extracted_coincident_element_ids(stk::CommSparse &commSparse,
                                           const std::map<int, stk::mesh::EntityIdVector> &extractedEntityIdsByProc);
void unpack_extracted_coincident_element_ids_and_delete_connecting_graph_edges(stk::CommSparse &commSparse, int procId, stk::mesh::Graph &graph);
void remove_edges_to_extracted_coincident_elements_on_other_procs(const std::map<int, stk::mesh::EntityIdVector> &extractedEntityIdsByProc,
                                                                  stk::mesh::Graph &graph,
                                                                  MPI_Comm comm);
}
}} // end namespaces stk mesh

#endif
