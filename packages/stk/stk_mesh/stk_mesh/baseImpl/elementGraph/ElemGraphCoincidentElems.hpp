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

class CoincidenceDetector
{
public:
    virtual ~CoincidenceDetector() {}
    virtual bool are_graph_edge_elements_coincident(const stk::mesh::GraphEdge &graphEdge) const = 0;
    virtual void report_coincident_sides(std::ostream &stream,
                                         const GraphEdgeVector& partiallyCoincidentSides) const {}
};

class CoincidentSideExtractor
{
public:
    CoincidentSideExtractor(stk::mesh::Graph &graph,
                            const std::vector<stk::topology> &topologies,
                            const CoincidenceDetector &detector)
    : m_graph(graph),
      m_topologies(topologies),
      m_detector(detector) {}

    SparseGraph extract_coincident_sides();
    void append_extracted_coincident_sides(const std::vector<impl::LocalId> &elemIds,
                                           stk::mesh::impl::SparseGraph &coincidentEdges);
private:
    CoincidentSideExtractor();

    void extract_coincident_sides(SparseGraph& extractedCoincidentSides, const CoincidenceDetector &detector);
    void extract_coincident_sides_for_element(stk::mesh::impl::LocalId elemId, GraphEdgeVector &partiallyCoincidentSides, const CoincidenceDetector &detector);
    void extract_coincident_sides_for_element(stk::mesh::impl::LocalId elemId, SparseGraph& extractedCoincidentSides, const CoincidenceDetector &detector);
    void delete_edges(const GraphEdgeVector& edgesToDelete);
    void add_edges(const GraphEdgeVector& edgesToDelete, SparseGraph& extractedCoincidentSides);

    stk::mesh::Graph &m_graph;
    const std::vector<stk::topology> &m_topologies;
    const CoincidenceDetector &m_detector;
};

void make_chosen_ids_in_parinfo_consistent_for_edges_with_coincident_elements(const stk::mesh::Graph &graph,
                                            stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                            const stk::mesh::impl::SparseGraph &extractedCoincidentElements,
                                            const IdMapper &idMapper,
                                            MPI_Comm comm);
}}} // end namespaces stk mesh

#endif
