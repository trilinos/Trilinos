#include "gtest/gtest.h"
#include <vector>
#include <mpi.h>
#include <stk_mesh/base/ElemElemGraph.hpp>

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
    EXPECT_EQ(1u, extractedCoincidentElements.size());
    for(int side = 0; side < elemDesc.numSides; side++)
        EXPECT_EQ(stk::mesh::GraphEdge(elemDesc.elem1, side, elemDesc.elem2, side), extractedCoincidentElements[elemDesc.elem1][side]);
}

void test_extracting_coincident_hex8s(stk::mesh::Graph &graph, const std::vector<stk::mesh::impl::CoincidentElementDescription> &elemDescs)
{
    std::vector<stk::topology> topologies = {stk::topology::HEX_8, stk::topology::HEX_8};
    stk::mesh::impl::SparseGraph extractedCoincidentElements = stk::mesh::impl::extract_coincident_sides(graph, topologies);

    expect_num_edges_remaining_per_element(graph, {0, 0});
    expect_coincident_elements_edges_were_extracted(extractedCoincidentElements, elemDescs[0]);
    EXPECT_TRUE(extractedCoincidentElements.find(1) == extractedCoincidentElements.end());
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
    expect_num_edges_remaining_per_element(graph, {1, 0, 1});
    expect_graph_edge_validity_per_element(graph, {true, false, true});
    expect_coincident_elements_edges_were_extracted(extractedCoincidentElements, {6, 0, 1});
    EXPECT_TRUE(extractedCoincidentElements.find(1) == extractedCoincidentElements.end());
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

void add_remote_edge_to_graph(stk::mesh::Graph &graph, stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges, const stk::mesh::GraphEdge &graphEdge, int otherProc)
{
    graph.add_edge(graphEdge);
    stk::mesh::impl::parallel_info parInfo(otherProc, 0, 1, stk::topology::HEX_8, true);
    parallelInfoForGraphEdges.insert_parallel_info_for_graph_edge(graphEdge, parInfo);
}

void make_graph_of_coincident_hex8s_with_adjacent_hex_in_parallel(stk::mesh::Graph &graph, stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges, MPI_Comm comm)
{
    if(stk::parallel_machine_rank(comm) == 0)
    {
        std::vector<stk::mesh::impl::CoincidentElementDescription> elemDescs = {{6, 0, 1}, {6, 1, 0}};
        make_graph_for_coincident_elements(graph, 2, elemDescs);
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(1, 3, -3, 4), 1);
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(0, 3, -3, 4), 1);
    }
    else
    {
        graph.add_new_element();
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(0, 4, -1, 3), 0);
        add_remote_edge_to_graph(graph, parallelInfoForGraphEdges, stk::mesh::GraphEdge(0, 4, -2, 3), 0);
    }
}

bool running_on_num_procs(MPI_Comm comm, int targetProcs)
{
    return stk::parallel_machine_size(comm) == targetProcs;
}

TEST(CoincidentElements, ExtractCoincidentHex8sWithAdjacentHexInParallel)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    if(running_on_num_procs(comm, 2))
    {
        stk::mesh::Graph graph;
        stk::mesh::ParallelInfoForGraphEdges parallelInfoForGraphEdges(stk::parallel_machine_rank(comm));
        make_graph_of_coincident_hex8s_with_adjacent_hex_in_parallel(graph, parallelInfoForGraphEdges, comm);

        std::vector<stk::topology> topologies = {stk::topology::HEX_8, stk::topology::HEX_8};
        stk::mesh::impl::SparseGraph extractedCoincidentElements = stk::mesh::impl::extract_coincident_sides(graph, topologies);

        std::map<int, std::vector<stk::mesh::impl::LocalId>> extractedIdsByProc =
                stk::mesh::impl::get_extracted_coincident_local_ids(graph, parallelInfoForGraphEdges, extractedCoincidentElements);
        if(stk::parallel_machine_rank(comm) == 0)
        {
            ASSERT_EQ(1u, extractedIdsByProc.size());
            EXPECT_EQ(1u, extractedIdsByProc[1].size());
            EXPECT_EQ(1, extractedIdsByProc[1][0]);
        }
        else
        {
            ASSERT_TRUE(extractedIdsByProc.empty());
        }

        std::map<int, stk::mesh::EntityIdVector> extractedEntityIdsByProc;
        if(stk::parallel_machine_rank(comm) == 0)
            extractedEntityIdsByProc[1] = {2};

        stk::mesh::impl::remove_edges_to_extracted_coincident_elements_on_other_procs(extractedEntityIdsByProc, graph, comm);

        if(stk::parallel_machine_rank(comm) == 0)
        {
            EXPECT_EQ(1u, graph.get_num_edges_for_element(0));
            EXPECT_EQ(0u, graph.get_num_edges_for_element(1));
        }
        else
        {
            EXPECT_EQ(1u, graph.get_num_edges_for_element(0));
        }
    }
}
}
