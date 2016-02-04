#ifndef STK_ELEM_ELEM_GRAPH_COINCIDENT_ELEMS_HPP
#define STK_ELEM_ELEM_GRAPH_COINCIDENT_ELEMS_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class Graph; } }
namespace stk { namespace mesh { class ParallelInfoForGraphEdges; } }

namespace stk
{
namespace mesh
{
namespace impl
{

struct CoincidentElementDescription
{
    int numSides;
    stk::mesh::impl::LocalId elem1;
    stk::mesh::impl::LocalId elem2;
};

typedef std::map<stk::mesh::impl::LocalId, std::vector<stk::mesh::GraphEdge>> SparseGraph;

class IdMapper
{
public:
    virtual ~IdMapper() { }
    virtual stk::mesh::EntityId localToGlobal(stk::mesh::impl::LocalId local) const = 0;
    virtual stk::mesh::impl::LocalId globalToLocal(stk::mesh::EntityId global) const = 0;
};

class CoincidentSideExtractor
{
public:
    CoincidentSideExtractor(stk::mesh::BulkData &bulkData,
                            stk::mesh::Graph &graph,
                            const std::vector<stk::topology> &topologies,
                            const stk::mesh::EntityVector &localIdToElementEntity,
                            const ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
    : m_bulkData(bulkData),
      m_graph(graph),
      m_topologies(topologies),
      m_localIdToElementEntity(localIdToElementEntity),
      m_parallelInfoForGraphEdges(parallelInfoForGraphEdges) {}

    SparseGraph extract_coincident_sides();
    void append_extracted_coincident_sides(const std::vector<impl::LocalId> &elemIds,
                                           stk::mesh::impl::SparseGraph &coincidentEdges);
private:
    CoincidentSideExtractor();

    void extract_partially_coincident_sides(SparseGraph& extractedCoincidentSides);
    stk::mesh::EntityVector get_side_nodes(const impl::LocalId elemId, const int side);
    bool are_graph_edge_elements_partially_coincident(const stk::mesh::GraphEdge &graphEdge);
    bool are_local_graph_edge_elements_partially_coincident(const stk::mesh::GraphEdge &graphEdge);
    bool are_parallel_graph_edge_elements_partially_coincident(const stk::mesh::GraphEdge &graphEdge);
    bool do_side_nodes_for_graph_edge_have_same_polarity(const stk::mesh::GraphEdge &graphEdge,
                                                         const stk::mesh::EntityVector &sideNodesElement1,
                                                         const stk::mesh::EntityVector &sideNodesElement2);
    void fill_partially_coincident_sides(stk::mesh::impl::LocalId elemId, GraphEdgeVector &partiallyCoincidentSides);
    int count_shared_sides(stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2);
    bool are_elements_coincident(int numSides, stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2);
    void fill_coincident_edges(size_t elemId, GraphEdgeVector &edgesToDelete);
    void delete_edges(const GraphEdgeVector& edgesToDelete);
    void add_edges(const GraphEdgeVector& edgesToDelete, SparseGraph& extractedCoincidentSides);
    void extract_coincident_sides_for_element(size_t elemId, SparseGraph& extractedCoincidentSides);

    stk::mesh::BulkData &m_bulkData;
    stk::mesh::Graph &m_graph;
    const std::vector<stk::topology> &m_topologies;
    const stk::mesh::EntityVector &m_localIdToElementEntity;
    const ParallelInfoForGraphEdges &m_parallelInfoForGraphEdges;
};

void make_chosen_ids_in_parinfo_consistent_for_edges_with_coincident_elements(const stk::mesh::Graph &graph,
                                            stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                            const stk::mesh::impl::SparseGraph &extractedCoincidentElements,
                                            const IdMapper &idMapper,
                                            MPI_Comm comm);
}}} // end namespaces stk mesh

#endif
