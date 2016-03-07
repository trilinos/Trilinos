#include <mpi.h>                        // for MPI_Comm, etc
#include <stddef.h>                     // for size_t
#include <map>                          // for map, map<>::mapped_type, etc
#include <ostream>                      // for basic_ostream::operator<<
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for EXPECT_EQ, TEST, etc
#include "stk_mesh/base/Types.hpp"      // for EntityIdVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc
#include <stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemGraphCoincidentElems.hpp>
#include "stk_mesh/baseImpl/elementGraph//ElemElemGraphImpl.hpp"

namespace
{

struct CoincidentElementDescription
{
    int numSides;
    stk::mesh::impl::LocalId elem1;
    stk::mesh::impl::LocalId elem2;
};

void add_remote_edge_to_graph(stk::mesh::Graph &graph,
                              stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                              const stk::mesh::GraphEdge &graphEdge,
                              int otherProc,
                              stk::mesh::EntityId chosenFaceId)
{
    graph.add_edge(graphEdge);
    stk::mesh::impl::ParallelInfo parInfo(otherProc, 0, chosenFaceId, stk::topology::HEX_8, true);
    parallelInfoForGraphEdges.insert_parallel_info_for_graph_edge(graphEdge, parInfo);
}

void make_graph_of_coincident_hex8s_with_adjacent_hex_in_parallel(stk::mesh::Graph &graph,
                                                                  stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                                                  MPI_Comm comm)
{
    if(stk::parallel_machine_rank(comm) == 0)
    {
        std::vector<CoincidentElementDescription> elemDescs = {{6, 0, 1}, {6, 1, 0}};
        graph.set_num_local_elements(2);
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(0, 3, -3, 4), 1, 13);
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(1, 3, -3, 4), 1, 23);
    }
    else
    {
        graph.add_new_element();
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(0, 4, -1, 3), 0, 13);
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(0, 4, -2, 3), 0, 23);
    }
}

bool running_on_num_procs(MPI_Comm comm, int targetProcs)
{
    return stk::parallel_machine_size(comm) == targetProcs;
}

class MockIdMapper : public stk::mesh::impl::IdMapper
{
public:
    MockIdMapper(const std::vector<stk::mesh::EntityId> &g) : globalIds(g) {}
    virtual ~MockIdMapper() {}
    virtual stk::mesh::EntityId localToGlobal(stk::mesh::impl::LocalId local) const
    {
        if(local<0)
            return -local;
        else
            return globalIds[local];
    }
    virtual stk::mesh::impl::LocalId globalToLocal(stk::mesh::EntityId global) const
    {
        auto iter = std::lower_bound(globalIds.begin(), globalIds.end(), global);
        if(iter == globalIds.end())
            return -1;
        return iter - globalIds.begin();
    }
private:
    std::vector<stk::mesh::EntityId> globalIds;
};

void expect_chosen_side_id_for_graph_edge(const stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                          const stk::mesh::GraphEdge &graphEdge,
                                          stk::mesh::EntityId expectedChosenSideId)
{
    const stk::mesh::impl::ParallelInfo &parInfo =
            parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
    EXPECT_EQ(expectedChosenSideId, parInfo.m_chosen_side_id);
}

void add_edges_between_elems(stk::mesh::impl::SparseGraph &coincidentElements, stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2, unsigned numSides)
{
    for(unsigned side=0; side<numSides; side++)
        coincidentElements.add_edge(stk::mesh::GraphEdge(elem1,side,elem2,side));
}

void add_coincident_edges_between_elems(stk::mesh::impl::SparseGraph &coincidentElements, stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2, unsigned numSides)
{
    add_edges_between_elems(coincidentElements, elem1, elem2, numSides);
    add_edges_between_elems(coincidentElements, elem2, elem1, numSides);
}

TEST(CoincidentElements, CorrectFaceId)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    if(running_on_num_procs(comm, 2))
    {
        stk::mesh::Graph graph;
        stk::mesh::ParallelInfoForGraphEdges parallelInfoForGraphEdges(stk::parallel_machine_rank(comm));
        make_graph_of_coincident_hex8s_with_adjacent_hex_in_parallel(graph, parallelInfoForGraphEdges, comm);

        const unsigned numSidesOfHex = 6;

        stk::mesh::impl::SparseGraph coincidentElements;
        if(stk::parallel_machine_rank(comm) == 0)
            add_coincident_edges_between_elems(coincidentElements, 0, 1, numSidesOfHex);

        std::vector<stk::mesh::EntityId> globalIds = {1, 2};
        if(stk::parallel_machine_rank(comm) == 1)
            globalIds = {3};

        MockIdMapper idMapper(globalIds);
        make_chosen_ids_in_parinfo_consistent_for_edges_with_coincident_elements(graph, parallelInfoForGraphEdges, coincidentElements, idMapper, comm);

        if(stk::parallel_machine_rank(comm) == 0)
        {
            stk::mesh::GraphEdge graphEdgeHex1Hex3(0, 3, -3, 4);
            expect_chosen_side_id_for_graph_edge(parallelInfoForGraphEdges, graphEdgeHex1Hex3, 13);
            stk::mesh::GraphEdge graphEdgeHex2Hex3(1, 3, -3, 4);
            expect_chosen_side_id_for_graph_edge(parallelInfoForGraphEdges, graphEdgeHex2Hex3, 13);
        }
        else
        {
            stk::mesh::GraphEdge graphEdgeHex3Hex1(0, 4, -1, 3);
            expect_chosen_side_id_for_graph_edge(parallelInfoForGraphEdges, graphEdgeHex3Hex1, 13);
            stk::mesh::GraphEdge graphEdgeHex3Hex2(0, 4, -2, 3);
            expect_chosen_side_id_for_graph_edge(parallelInfoForGraphEdges, graphEdgeHex3Hex2, 13);
        }
    }
}

}
