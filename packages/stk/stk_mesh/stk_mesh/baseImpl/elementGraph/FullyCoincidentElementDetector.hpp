#ifndef FULLY_COINCIDENT_ELEMENT_DETECTOR_HPP
#define FULLY_COINCIDENT_ELEMENT_DETECTOR_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include <stk_topology/topology.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

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
        return are_elements_fully_coincident(static_cast<int>(m_topologies[graphEdge.elem1].num_sides()),
                                             graphEdge.elem1, graphEdge.elem2);
    }

    virtual void report_coincident_sides(std::ostream &stream, const GraphEdgeVector& coincidentSides) const
    {
        if(coincidentSides.size() > 0)
        {
            std::ostringstream os;
            os << "There are " << coincidentSides.size() << " co-incident edges" << std::endl;
            for(const auto &graphEdge : coincidentSides)
            {
                os << "     (" << graphEdge.elem1
                        << "," << graphEdge.side1
                        << ")  is co-incident with "
                        << "(" << graphEdge.elem2
                        << "," << graphEdge.side2
                        << ")"
                        << std::endl;
            }
            stream << os.str();
        }
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


}
}
}

#endif
