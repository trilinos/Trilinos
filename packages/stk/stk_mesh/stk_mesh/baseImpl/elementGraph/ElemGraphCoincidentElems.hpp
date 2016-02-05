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
    virtual void report_partially_coincident_sides(std::ostream &stream,
                                                   const GraphEdgeVector& partiallyCoincidentSides) const {}
};

class FullyCoincidentElementDetector: public CoincidenceDetector
{
public:
    FullyCoincidentElementDetector(const stk::mesh::Graph &graph,
                                    const std::vector<stk::topology> &topologies)
    : m_graph(graph),
      m_topologies(topologies) {}

    virtual ~FullyCoincidentElementDetector() {}
    virtual bool are_graph_edge_elements_coincident(const stk::mesh::GraphEdge &graphEdge) const
    {
        return are_elements_fully_coincident(static_cast<int>(m_topologies[graphEdge.elem1].num_sides()), graphEdge.elem1, graphEdge.elem2);
    }
private:

    int count_shared_sides(stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2) const
    {
        int numSharedSides = 0;
        for(size_t i=0; i < m_graph.get_num_edges_for_element(elem1); i++)
        {
            const stk::mesh::GraphEdge &graphEdge = m_graph.get_edge_for_element(elem1, i);
            if(graphEdge.elem2 == elem2)
                numSharedSides++;
        }
        return numSharedSides;
    }

    bool are_elements_fully_coincident(int numSides, stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2) const
    {
        int numSharedSides = count_shared_sides(elem1, elem2);
        return (numSharedSides == numSides);
    }

    const stk::mesh::Graph &m_graph;
    const std::vector<stk::topology> &m_topologies;
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

    void extract_partially_coincident_sides(SparseGraph& extractedCoincidentSides);
    void fill_partially_coincident_sides(stk::mesh::impl::LocalId elemId, GraphEdgeVector &partiallyCoincidentSides);
    int count_shared_sides(stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2);
    bool are_elements_fully_coincident(int numSides, stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2);
    void fill_fully_coincident_edges(size_t elemId, GraphEdgeVector &edgesToDelete);
    void delete_edges(const GraphEdgeVector& edgesToDelete);
    void add_edges(const GraphEdgeVector& edgesToDelete, SparseGraph& extractedCoincidentSides);
    void extract_fully_coincident_sides_for_element(size_t elemId, SparseGraph& extractedCoincidentSides);

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
