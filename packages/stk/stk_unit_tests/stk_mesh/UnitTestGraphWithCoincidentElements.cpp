#include "gtest/gtest.h"
#include <vector>
#include <mpi.h>
#include <stk_mesh/base/GraphEdgeData.hpp>
#include <stk_mesh/base/ElemGraphCoincidentElems.hpp>

namespace
{
void make_edges_for_coincident_elems(stk::mesh::Graph &graph, const stk::mesh::impl::CoincidentElementDescription &elemDesc)
{
    for(int side=0; side < elemDesc.numSides; side++)
        graph.add_edge(stk::mesh::GraphEdge(elemDesc.elem1, side, elemDesc.elem2, side));
}

void make_graph_for_coincident_elements(stk::mesh::Graph &graph, size_t numLocalElems, const std::vector<stk::mesh::impl::CoincidentElementDescription> &elemDescs)
{
    graph.set_num_local_elements(numLocalElems);
    for(const stk::mesh::impl::CoincidentElementDescription &elemDesc : elemDescs)
        make_edges_for_coincident_elems(graph, elemDesc);
}

void add_symmetric_edges(stk::mesh::Graph &graph, const stk::mesh::GraphEdge &graphEdge)
{
    graph.add_edge(graphEdge);
    graph.add_edge(create_symmetric_edge(graphEdge));
}

void add_adjacent_hex_to_graph(stk::mesh::Graph &graph, int sideFromStacked, int sideFromNew)
{
    graph.add_new_element();
    add_symmetric_edges(graph, stk::mesh::GraphEdge(0, sideFromStacked, 2, sideFromNew));
    add_symmetric_edges(graph, stk::mesh::GraphEdge(1, sideFromStacked, 2, sideFromNew));
}

void make_graph_of_coincident_hex8s_with_adjacent_hex(stk::mesh::Graph &graph, int sideFromStacked, int sideFromNew)
{
    std::vector<stk::mesh::impl::CoincidentElementDescription> elemDescs = {{6, 0, 1}, {6, 1, 0}};
    make_graph_for_coincident_elements(graph, 2, elemDescs);
    add_adjacent_hex_to_graph(graph, sideFromStacked, sideFromNew);
}

void expect_num_edges_remaining_per_element(const stk::mesh::Graph &graph, const std::vector<size_t> &goldNumEdgesRemaining)
{
    for(size_t elem = 0; elem < goldNumEdgesRemaining.size(); elem++)
        EXPECT_EQ(goldNumEdgesRemaining[elem], graph.get_num_edges_for_element(elem));
}

void expect_coincident_elements_edges_were_extracted(stk::mesh::impl::SparseGraph &extractedCoincidentElements,
                                                     const stk::mesh::impl::CoincidentElementDescription &elemDesc)
{
    EXPECT_EQ(2u, extractedCoincidentElements.size());
    for(int side = 0; side < elemDesc.numSides; side++)
    {
        EXPECT_EQ(stk::mesh::GraphEdge(elemDesc.elem1, side, elemDesc.elem2, side), extractedCoincidentElements[elemDesc.elem1][side]);
        EXPECT_EQ(stk::mesh::GraphEdge(elemDesc.elem2, side, elemDesc.elem1, side), extractedCoincidentElements[elemDesc.elem2][side]);
    }
}

void test_extracting_coincident_hex8s(stk::mesh::Graph &graph, const std::vector<stk::mesh::impl::CoincidentElementDescription> &elemDescs)
{
    std::vector<stk::topology> topologies = {stk::topology::HEX_8, stk::topology::HEX_8};
    stk::mesh::impl::SparseGraph extractedCoincidentElements = stk::mesh::impl::extract_coincident_sides(graph, topologies);

    expect_num_edges_remaining_per_element(graph, {0, 0});
    expect_coincident_elements_edges_were_extracted(extractedCoincidentElements, elemDescs[0]);
//    EXPECT_TRUE(extractedCoincidentElements.find(1) == extractedCoincidentElements.end());
}

TEST(CoincidentElements, ExtractCoincidentHex8s)
{
    stk::mesh::Graph graph;
    std::vector<stk::mesh::impl::CoincidentElementDescription> elemDescs = {{6, 0, 1}, {6, 1, 0}};
    make_graph_for_coincident_elements(graph, 2, elemDescs);
    test_extracting_coincident_hex8s(graph, elemDescs);
}

void expect_graph_edge_validity_per_element(const stk::mesh::Graph &graph, const std::vector<bool> &goldValidityPerElem)
{
    for(size_t elem = 0; elem < goldValidityPerElem.size(); elem++)
        EXPECT_EQ(goldValidityPerElem[elem], graph.is_valid(elem)) << elem;
}

void expect_coincident_edges_removed_others_remain(const stk::mesh::Graph &graph, stk::mesh::impl::SparseGraph &extractedCoincidentElements)
{
    expect_num_edges_remaining_per_element(graph, {1, 1, 2});
    expect_graph_edge_validity_per_element(graph, {true, true, true});
    expect_coincident_elements_edges_were_extracted(extractedCoincidentElements, {6, 0, 1});
//    EXPECT_TRUE(extractedCoincidentElements.find(1) == extractedCoincidentElements.end());
}

void test_extracting_coincident_hex8s_with_adjacent_hex(stk::mesh::Graph &graph)
{
    std::vector<stk::topology> topologies = {stk::topology::HEX_8, stk::topology::HEX_8, stk::topology::HEX_8};
    stk::mesh::impl::SparseGraph extractedCoincidentElements = stk::mesh::impl::extract_coincident_sides(graph, topologies);
    expect_coincident_edges_removed_others_remain(graph, extractedCoincidentElements);
}

TEST(CoincidentElements, ExtractCoincidentHex8sWithAdjacentHex)
{
    stk::mesh::Graph graph;
    int sideFromStacked = 4, sideFromNew = 3;
    make_graph_of_coincident_hex8s_with_adjacent_hex(graph, sideFromStacked, sideFromNew);
    test_extracting_coincident_hex8s_with_adjacent_hex(graph);
}

void add_remote_edge_to_graph(stk::mesh::Graph &graph,
                              stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                              const stk::mesh::GraphEdge &graphEdge,
                              int otherProc,
                              stk::mesh::EntityId chosenFaceId)
{
    graph.add_edge(graphEdge);
    stk::mesh::impl::parallel_info parInfo(otherProc, 0, chosenFaceId, stk::topology::HEX_8, true);
    parallelInfoForGraphEdges.insert_parallel_info_for_graph_edge(graphEdge, parInfo);
}

void make_graph_of_coincident_hex8s_with_adjacent_hex_in_parallel(stk::mesh::Graph &graph,
                                                                  stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                                                  MPI_Comm comm)
{
    if(stk::parallel_machine_rank(comm) == 0)
    {
        std::vector<stk::mesh::impl::CoincidentElementDescription> elemDescs = {{6, 0, 1}, {6, 1, 0}};
        make_graph_for_coincident_elements(graph, 2, elemDescs);
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
    const stk::mesh::impl::parallel_info &parInfo =
            parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
    EXPECT_EQ(expectedChosenSideId, parInfo.m_chosen_side_id);
}

TEST(CoincidentElements, CorrectFaceId)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    if(running_on_num_procs(comm, 2))
    {
        stk::mesh::Graph graph;
        stk::mesh::ParallelInfoForGraphEdges parallelInfoForGraphEdges(stk::parallel_machine_rank(comm));
        make_graph_of_coincident_hex8s_with_adjacent_hex_in_parallel(graph, parallelInfoForGraphEdges, comm);

        std::vector<stk::topology> topologies = {stk::topology::HEX_8, stk::topology::HEX_8};
        stk::mesh::impl::SparseGraph extractedCoincidentElements = stk::mesh::impl::extract_coincident_sides(graph, topologies);

        std::vector<stk::mesh::EntityId> globalIds = {1, 2};
        if(stk::parallel_machine_rank(comm) == 1)
            globalIds = {3};

        MockIdMapper idMapper(globalIds);
        choose_face_id_for_coincident_elements(graph, parallelInfoForGraphEdges, extractedCoincidentElements, idMapper, comm);

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
