#ifndef STK_ELEM_GRAPH_SHELL_CONNECTIONS_HPP
#define STK_ELEM_GRAPH_SHELL_CONNECTIONS_HPP

#include <vector>
#include <mpi.h>
#include <stk_mesh/base/GraphEdgeData.hpp>
#include <stk_mesh/base/ElemElemGraphImpl.hpp>

namespace stk
{
namespace mesh
{

struct GraphInfo
{
    GraphInfo(stk::mesh::Graph &g, stk::mesh::ParallelInfoForGraphEdges &p, std::vector<stk::topology> &e) :
            graph(g),
            parGraphInfo(p),
            elementTopologies(e)
    {
    }
    stk::mesh::Graph &graph;
    stk::mesh::ParallelInfoForGraphEdges &parGraphInfo;
    std::vector<stk::topology> &elementTopologies;
};

void remove_graph_edges_blocked_by_shell(GraphInfo &graphInfo);

}
}

#endif
