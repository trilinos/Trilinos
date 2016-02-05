#ifndef BULK_DATA_COINCIDENCE_DETECTOR_HPP
#define BULK_DATA_COINCIDENCE_DETECTOR_HPP

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

class BulkDataCoincidenceDetector : public CoincidenceDetector
{
public:
    BulkDataCoincidenceDetector(stk::mesh::BulkData &bulkData,
                                const std::vector<stk::topology> &topologies,
                                const stk::mesh::EntityVector &localIdToElementEntity,
                                const ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
        : m_bulkData(bulkData),
          m_topologies(topologies),
          m_localIdToElementEntity(localIdToElementEntity),
          m_parallelInfoForGraphEdges(parallelInfoForGraphEdges) {}

    virtual ~BulkDataCoincidenceDetector() {}
    virtual bool are_graph_edge_elements_coincident(const stk::mesh::GraphEdge &graphEdge) const
    {
        if(!impl::is_local_element(graphEdge.elem2))
            return are_parallel_graph_edge_elements_partially_coincident(graphEdge);

        return are_local_graph_edge_elements_partially_coincident(graphEdge);
    }

    virtual void report_coincident_sides(std::ostream &stream,
                                         const GraphEdgeVector& coincidentSides)
    {
#if 0
        if(coincidentSides.size() > 0)
        {
            std::ostringstream os;
            os << "There are " << coincidentSides.size() << " co-incident edges" << std::endl;
            for(const auto &graphEdge : coincidentSides)
            {
                os << "     (" << m_bulkData.identifier(m_localIdToElementEntity[graphEdge.elem1])
                       << "," << graphEdge.side1
                       << ")  is co-incident with "
                       << "(" << m_bulkData.identifier(m_localIdToElementEntity[graphEdge.elem2])
                       << "," << graphEdge.side2
                       << ")"
                       << std::endl;
            }
            stream << os.str();
        }
#else
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
#endif
    }

private:
    BulkDataCoincidenceDetector();

    stk::mesh::EntityVector get_side_nodes(const impl::LocalId elemId, const int side) const
    {
        unsigned numSideNodes = m_topologies[elemId].side_topology(side).num_nodes();
        stk::mesh::EntityVector sideNodes(numSideNodes);
        stk::mesh::Entity element = m_localIdToElementEntity[elemId];
        m_topologies[elemId].side_nodes(m_bulkData.begin_nodes(element), side, sideNodes.begin());
        return sideNodes;
    }

    bool are_parallel_graph_edge_elements_partially_coincident(const stk::mesh::GraphEdge &graphEdge) const
    {
        const impl::parallel_info& pInfo = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
        stk::topology remoteSideTopology = pInfo.m_remote_element_toplogy.side_topology(graphEdge.side2);
        bool same_polarity = static_cast<unsigned>(pInfo.m_permutation) < remoteSideTopology.num_positive_permutations();
        return false && same_polarity;
    }

    bool do_side_nodes_for_graph_edge_have_same_polarity(const stk::mesh::GraphEdge &graphEdge,
                                                         const stk::mesh::EntityVector &sideNodesElement1,
                                                         const stk::mesh::EntityVector &sideNodesElement2) const
    {
        stk::topology side1Topology = m_topologies[graphEdge.elem1].side_topology(graphEdge.side1);
        std::pair<bool,unsigned> result = side1Topology.equivalent(sideNodesElement1, sideNodesElement2);
        bool same_polarity = result.second < side1Topology.num_positive_permutations();
        return result.first && same_polarity;
    }

    bool are_local_graph_edge_elements_partially_coincident(const stk::mesh::GraphEdge &graphEdge) const
    {
        stk::mesh::EntityVector sideNodesElement1 = get_side_nodes(graphEdge.elem1, graphEdge.side1);
        stk::mesh::EntityVector sideNodesElement2 = get_side_nodes(graphEdge.elem2, graphEdge.side2);

        if(sideNodesElement1.size() != sideNodesElement2.size())
            return false;

        return do_side_nodes_for_graph_edge_have_same_polarity(graphEdge, sideNodesElement1, sideNodesElement2);
    }

    stk::mesh::BulkData &m_bulkData;
    const std::vector<stk::topology> &m_topologies;
    const stk::mesh::EntityVector &m_localIdToElementEntity;
    const ParallelInfoForGraphEdges &m_parallelInfoForGraphEdges;
};

}
}
}

#endif
