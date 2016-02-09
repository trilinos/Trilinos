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
        return stk::mesh::impl::are_graph_edge_elements_fully_coincident(m_graph, m_topologies, graphEdge);
    }

    virtual void report_coincident_sides(std::ostream &stream, const GraphEdgeVector& coincidentSides) const
    {
        if(coincidentSides.size() > 0)
        {
            std::ostringstream os;
            os << "FullyCoincidentElementDetector: there are " << coincidentSides.size() << " co-incident edges" << std::endl;
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
    const stk::mesh::Graph &m_graph;
    const std::vector<stk::topology> &m_topologies;
};


}
}
}

#endif
