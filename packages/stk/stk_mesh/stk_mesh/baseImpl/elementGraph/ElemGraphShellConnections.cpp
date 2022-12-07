#include "ElemGraphShellConnections.hpp"
#include <vector>
#include "mpi.h"
#include "GraphEdgeData.hpp"
#include "ElemElemGraphImpl.hpp"

namespace stk
{
namespace mesh
{

bool is_other_element_shell(const stk::mesh::GraphEdge& graphEdge,
                            const std::vector<stk::topology>& elementTopologies,
                            const stk::mesh::ParallelInfoForGraphEdges& parGraphInfo)
{
    stk::topology topo;
    if(stk::mesh::impl::is_local_element(graphEdge.elem2()))
        topo = elementTopologies[graphEdge.elem2()];
    else
        topo = parGraphInfo.get_parallel_info_for_graph_edge(graphEdge).m_remote_element_topology;
    return stk::mesh::impl::is_shell_or_beam2(topo);
}

bool is_this_element_shell(const stk::mesh::GraphEdge& graphEdge,
                          const std::vector<stk::topology>& elementTopologies)
{
    stk::topology topo = elementTopologies[graphEdge.elem1()];
    return stk::mesh::impl::is_shell_or_beam2(topo);
}

constexpr int MAX_NUM_SIDE_CONNECTIONS = 12;

class SideConnections
{
public:
    SideConnections(size_t numSides_);
    std::vector<int> get_sides_connected_to_shell_and_nonshell(const GraphInfo &graphInfo, size_t localId);

private:
    void record_side_connections_for_local_element(const GraphInfo &graphInfo, size_t localId);
    void record_side_connections_for_graph_edge(const stk::mesh::GraphEdge& graphEdge,
                                                const stk::mesh::ParallelInfoForGraphEdges& parGraphInfo,
                                                const std::vector<stk::topology>& elementTopologies);
    void fill_sides_connected_to_shell_and_nonshell(std::vector<int>& sidesConnectedToShellAndToNonShell);
    bool does_side_have_both_connection_to_shell_and_to_nonshell(int side) const;

private:
    size_t numSides;
    bool foundShellViaSide[MAX_NUM_SIDE_CONNECTIONS];
    bool foundSomethingUnShellLikeViaSide[MAX_NUM_SIDE_CONNECTIONS];
};

SideConnections::SideConnections(size_t numSides_) :
        numSides(numSides_)
{
  for(int i=0; i<MAX_NUM_SIDE_CONNECTIONS; ++i) {
    foundShellViaSide[i] = false;
    foundSomethingUnShellLikeViaSide[i] = false;
  }
}

std::vector<int> SideConnections::get_sides_connected_to_shell_and_nonshell(const GraphInfo &graphInfo, size_t localId)
{
    record_side_connections_for_local_element(graphInfo, localId);
    std::vector<int> sidesConnectedToShellAndToNonShell;
    fill_sides_connected_to_shell_and_nonshell(sidesConnectedToShellAndToNonShell);
    return sidesConnectedToShellAndToNonShell;
}

void SideConnections::record_side_connections_for_local_element(const GraphInfo &graphInfo, size_t localId)
{
    for(const stk::mesh::GraphEdge &graphEdge : graphInfo.graph.get_edges_for_element(localId))
        record_side_connections_for_graph_edge(graphEdge, graphInfo.parGraphInfo, graphInfo.elementTopologies);
}

void SideConnections::record_side_connections_for_graph_edge(const stk::mesh::GraphEdge& graphEdge,
                                            const stk::mesh::ParallelInfoForGraphEdges& parGraphInfo,
                                            const std::vector<stk::topology>& elementTopologies)
{
    if(is_other_element_shell(graphEdge, elementTopologies, parGraphInfo))
        foundShellViaSide[graphEdge.side1()] = true;
    else
        foundSomethingUnShellLikeViaSide[graphEdge.side1()] = true;
}

void SideConnections::fill_sides_connected_to_shell_and_nonshell(std::vector<int>& sidesConnectedToShellAndToNonShell)
{
    for(size_t side = 0; side < numSides; side++)
        if(does_side_have_both_connection_to_shell_and_to_nonshell(side))
            sidesConnectedToShellAndToNonShell.push_back(side);
}

bool SideConnections::does_side_have_both_connection_to_shell_and_to_nonshell(int side) const
{
    return (foundShellViaSide[side] && foundSomethingUnShellLikeViaSide[side]);
}


void fill_non_shell_graph_edges_to_delete(GraphInfo &graphInfo, const stk::mesh::impl::ElementSidePair &elementSidePair, std::vector<GraphEdge>& edgesToDelete)
{
    std::vector<GraphEdge> pllEdgesToDelete;
    for(const stk::mesh::GraphEdge& graphEdge : graphInfo.graph.get_edges_for_element_side(elementSidePair.first, elementSidePair.second))
    {
        bool thisIsShell = is_this_element_shell(graphEdge, graphInfo.elementTopologies);
        bool thatIsShell = is_other_element_shell(graphEdge, graphInfo.elementTopologies, graphInfo.parGraphInfo);

        if(!thisIsShell && !thatIsShell)
        {
            if(!impl::is_local_element(graphEdge.elem2())) {
                pllEdgesToDelete.push_back(graphEdge);
            }
            edgesToDelete.push_back(graphEdge);
        }
    }
    graphInfo.parGraphInfo.erase_edges(pllEdgesToDelete);
}

void remove_graph_edges_blocked_by_shell(GraphInfo &graphInfo)
{
    std::vector<GraphEdge> edgesToDelete;
    for(size_t localId = 0; localId < graphInfo.graph.get_num_elements_in_graph(); localId++)
    {
        SideConnections sideConnectionsForElement(graphInfo.elementTopologies[localId].num_sides());
        for(int side : sideConnectionsForElement.get_sides_connected_to_shell_and_nonshell(graphInfo, localId))
            fill_non_shell_graph_edges_to_delete(graphInfo, stk::mesh::impl::ElementSidePair(localId, side), edgesToDelete);

        if (edgesToDelete.size() > 0)
        {
            std::sort(edgesToDelete.begin(), edgesToDelete.end(), GraphEdgeLessByElem1());
            graphInfo.graph.delete_sorted_edges(edgesToDelete);
            edgesToDelete.clear();
        }
    }
}

}
}
