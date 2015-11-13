#include "gtest/gtest.h"
#include <mpi.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_util/parallel/CommSparse.hpp>

namespace
{

bool running_on_num_procs(MPI_Comm comm, int targetProcs)
{
    return stk::parallel_machine_size(comm) == targetProcs;
}

int get_other_proc(MPI_Comm comm)
{
    return stk::parallel_machine_size(comm) - stk::parallel_machine_rank(comm) - 1;
}

class CoincidentElements: public stk::unit_test_util::MeshFixture
{
protected:
    void declare_num_coincident_elements(unsigned numElemsToCreate, const stk::mesh::EntityIdVector &nodes, stk::mesh::Part &part)
    {
        get_bulk().modification_begin();
        for(unsigned i = 0; i < numElemsToCreate; i++)
            stk::mesh::declare_element(get_bulk(), part, i+1, nodes);
        get_bulk().modification_end();
    }
    void declare_num_coincident_elements_round_robin(unsigned numElemsToCreate, const stk::mesh::EntityIdVector &nodes, stk::mesh::Part &part)
    {
        get_bulk().modification_begin();
        for(unsigned i = 0; i < numElemsToCreate; i++)
            if(static_cast<int>(i) % get_bulk().parallel_size() == get_bulk().parallel_rank())
                stk::mesh::declare_element(get_bulk(), part, i+1, nodes);
        for(const stk::mesh::EntityId nodeId : nodes)
            get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, nodeId), get_other_proc(get_comm()));
        get_bulk().modification_end();
    }
    void run_skin_mesh()
    {
        stk::mesh::ElemElemGraph elementGraph(get_bulk(), get_meta().universal_part());
        elementGraph.skin_mesh({});
    }
    void expect_faces_connected_to_num_elements_locally(const std::vector<unsigned> &goldNumConnectedElems)
    {
        stk::mesh::EntityVector faces;
        stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
        ASSERT_EQ(goldNumConnectedElems.size(), faces.size());
        for(size_t i=0; i<faces.size(); i++)
            EXPECT_EQ(goldNumConnectedElems[i], get_bulk().num_elements(faces[i])) << i;
    }
};

class CoincidentQuad4Shells : public CoincidentElements
{
protected:
    void skin_num_coincident_shells(unsigned numElemsToCreate)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        create_stacked_shells(numElemsToCreate);
        run_skin_mesh();
    }
    void create_stacked_shells(unsigned numElemsToCreate)
    {
        stk::mesh::Part &shellPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUADRILATERAL_4);
        stk::mesh::EntityIdVector nodes = {1, 2, 3, 4};
        declare_num_coincident_elements_round_robin(numElemsToCreate, nodes, shellPart);
    }

    void expect_local_face_shared_to_other_proc(const std::vector<int> &goldSharedProcs, const stk::mesh::EntityIdVector &goldFaceIds)
    {
        stk::mesh::EntityVector faces;
        stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
        ASSERT_EQ(goldSharedProcs.size(), faces.size());
        for(size_t i=0; i<faces.size(); i++)
        {
            EXPECT_EQ(goldFaceIds[i], get_bulk().identifier(faces[i]));

            std::vector<int> procs;
            get_bulk().comm_shared_procs(get_bulk().entity_key(faces[i]), procs);
            ASSERT_EQ(1u, procs.size());
            EXPECT_EQ(goldSharedProcs[i], procs[0]);
        }
    }
};

TEST_F(CoincidentQuad4Shells, Skin2)
{
    if(running_on_num_procs(get_comm(), 1))
    {
        unsigned numElemsToCreate = 2;
        skin_num_coincident_shells(numElemsToCreate);
        expect_faces_connected_to_num_elements_locally({numElemsToCreate, numElemsToCreate});
    }
}

TEST_F(CoincidentQuad4Shells, Skin2InParallel)
{
    if(running_on_num_procs(get_comm(), 2))
    {
        skin_num_coincident_shells(2);
        expect_faces_connected_to_num_elements_locally({1, 1});
        int otherProc = get_other_proc(get_comm());
        expect_local_face_shared_to_other_proc({otherProc, otherProc}, {1, 2});
    }
}

TEST_F(CoincidentQuad4Shells, Skin3)
{
    if(running_on_num_procs(get_comm(), 1))
    {
        unsigned numElemsToCreate = 3;
        skin_num_coincident_shells(numElemsToCreate);
        expect_faces_connected_to_num_elements_locally({numElemsToCreate, numElemsToCreate});
    }
}

class CoincidentHex8s : public CoincidentElements
{
protected:
    void create_two_coincident_hexes()
    {
        block1 = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
        stk::mesh::EntityIdVector nodes = {1, 2, 3, 4, 5, 6, 7, 8};
        declare_num_coincident_elements(2, nodes, *block1);
    }
    stk::mesh::Part *block1;
};
TEST_F(CoincidentHex8s, SkinMesh)
{
    if(running_on_num_procs(get_comm(), 1))
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        create_two_coincident_hexes();
        run_skin_mesh();
        expect_faces_connected_to_num_elements_locally({2, 2, 2, 2, 2, 2});
    }
}

class CoincidentHex8sWithAdjacentHex : public CoincidentHex8s
{
protected:
    void create_coincident_hex8s_with_adjacent_hex()
    {
        block2 = &get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        create_two_coincident_hexes();
        create_adjacent_hex();
    }
    void create_adjacent_hex()
    {
        stk::mesh::EntityIdVector nodes = {5, 6, 7, 8, 9, 10, 11, 12};
        get_bulk().modification_begin();
        stk::mesh::declare_element(get_bulk(), *block2, 3, nodes);
        get_bulk().modification_end();
    }
    stk::mesh::Part *block2;
};
TEST_F(CoincidentHex8sWithAdjacentHex, SkinMesh)
{
    if(running_on_num_procs(get_comm(), 1))
    {
        create_coincident_hex8s_with_adjacent_hex();

        stk::mesh::ElemElemGraph elementGraph(get_bulk(), get_meta().universal_part());
        elementGraph.skin_mesh({});

        expect_faces_connected_to_num_elements_locally({1, 1, 1, 1, 1, 2, 2, 2, 2, 2});
    }
}
TEST_F(CoincidentHex8sWithAdjacentHex, SkinMeshWithAirSelectorOnAdjacentHex)
{
    if(running_on_num_procs(get_comm(), 1))
    {
        create_coincident_hex8s_with_adjacent_hex();

        stk::mesh::Selector air = *block2;
        stk::mesh::ElemElemGraph elementGraph(get_bulk(), *block1, &air);
        elementGraph.skin_mesh( {});

        expect_faces_connected_to_num_elements_locally({2, 2, 2, 2, 2, 3});
    }
}
TEST_F(CoincidentHex8sWithAdjacentHex, SkinMeshWithAirSelectorOnCoincidentHexes)
{
    if(running_on_num_procs(get_comm(), 1))
    {
        create_coincident_hex8s_with_adjacent_hex();

        stk::mesh::Selector air = *block1;
        stk::mesh::ElemElemGraph elementGraph(get_bulk(), *block2, &air);
        elementGraph.skin_mesh( {});

        expect_faces_connected_to_num_elements_locally({1, 1, 1, 1, 1, 3});
    }
}

class CoincidentHex8sWithAdjacentHexInParallel : public CoincidentHex8sWithAdjacentHex
{
protected:
    void create_coincident_hex8s_with_adjacent_hex_on_2_procs()
    {
        block1 = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
        block2 = &get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        get_bulk().modification_begin();
        if(get_bulk().parallel_rank() == 0)
        {
            stk::mesh::EntityIdVector nodes = {1, 2, 3, 4, 5, 6, 7, 8};
            stk::mesh::declare_element(get_bulk(), *block1, 1, nodes);
            stk::mesh::declare_element(get_bulk(), *block1, 2, nodes);
        }
        else
            stk::mesh::declare_element(get_bulk(), *block2, 3, {5, 6, 7, 8, 9, 10, 11, 12});

        for(const stk::mesh::EntityId nodeId : {5, 6, 7, 8})
            get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, nodeId), get_other_proc(get_comm()));
        get_bulk().modification_end();
    }
    stk::mesh::Part *block2;
};
TEST_F(CoincidentHex8sWithAdjacentHexInParallel, SkinMeshWithAirSelectorOnAdjacentHex)
{
    if(running_on_num_procs(get_comm(), 2))
    {
        create_coincident_hex8s_with_adjacent_hex_on_2_procs();

        stk::mesh::Selector air = *block2;
        stk::mesh::ElemElemGraph elementGraph(get_bulk(), *block1, &air);
        elementGraph.skin_mesh({});

        if(get_bulk().parallel_rank() == 0)
            expect_faces_connected_to_num_elements_locally({2, 2, 2, 2, 2, 2});
        else
            expect_faces_connected_to_num_elements_locally({1});
    }
}
TEST_F(CoincidentHex8sWithAdjacentHexInParallel, SkinMeshWithAirSelectorOnCoincidentHexes)
{
    if(running_on_num_procs(get_comm(), 2))
    {
        create_coincident_hex8s_with_adjacent_hex_on_2_procs();

        stk::mesh::Selector air = *block1;
        stk::mesh::ElemElemGraph elementGraph(get_bulk(), *block2, &air);
        elementGraph.skin_mesh({});

        if(get_bulk().parallel_rank() == 0)
            expect_faces_connected_to_num_elements_locally({2});
        else
            expect_faces_connected_to_num_elements_locally({1, 1, 1, 1, 1, 1});
    }
}

void compare_exposed_sides(const std::vector<stk::mesh::CoincidentElementConnection> &goldExposedSides,
                           const std::vector<stk::mesh::CoincidentElementConnection> &exposedCoincidentSides)
{
    ASSERT_EQ(goldExposedSides.size(), exposedCoincidentSides.size());
    for(size_t i=0; i<exposedCoincidentSides.size();i++)
    {
        EXPECT_EQ(goldExposedSides[i], exposedCoincidentSides[i]);
    }
}

void expect_all_sides_exposed(const stk::mesh::Graph &graph, size_t numSides, stk::mesh::impl::LocalId elem1, std::vector<stk::mesh::impl::LocalId> connectedElems)
{
    std::vector<stk::mesh::CoincidentElementConnection> goldCoincidentElementConnections;
    for(size_t side=0; side < numSides; side++)
        for(stk::mesh::impl::LocalId connectedElem : connectedElems)
            goldCoincidentElementConnections.push_back(stk::mesh::CoincidentElementConnection(elem1, side, connectedElem, side));

    std::vector<stk::mesh::CoincidentElementConnection> exposedCoincidentSides = stk::mesh::impl::get_exposed_coincident_sides(graph, numSides, elem1);
    compare_exposed_sides(goldCoincidentElementConnections, exposedCoincidentSides);
}

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

void test_graph_for_coincident_elements(stk::mesh::Graph &graph, const std::vector<stk::mesh::impl::CoincidentElementDescription> &elemDescs)
{
    for(const stk::mesh::impl::CoincidentElementDescription &elemDesc : elemDescs)
        expect_all_sides_exposed(graph, elemDesc.numSides, elemDesc.elem1, {elemDesc.elem2});
}

TEST(CoincidentElements, DetectCoincidentHex8s)
{
    stk::mesh::Graph graph;
    std::vector<stk::mesh::impl::CoincidentElementDescription> elemDescs = {{6, 0, 1}, {6, 1, 0}};
    make_graph_for_coincident_elements(graph, 2, elemDescs);
    test_graph_for_coincident_elements(graph, elemDescs);
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

TEST(CoincidentElements, CoincidentHex8sWithAdjacentHex)
{
    stk::mesh::Graph graph;
    int sideFromStacked = 4, sideFromNew = 3;
    make_graph_of_coincident_hex8s_with_adjacent_hex(graph, sideFromStacked, sideFromNew);

    size_t numSidesHex8 = 6;
    std::vector<stk::mesh::CoincidentElementConnection> exposedCoincidentSides = stk::mesh::impl::get_exposed_coincident_sides(graph, numSidesHex8, 0);
    ASSERT_EQ(5u, exposedCoincidentSides.size());
    exposedCoincidentSides = stk::mesh::impl::get_exposed_coincident_sides(graph, numSidesHex8, 1);
    ASSERT_EQ(5u, exposedCoincidentSides.size());
    exposedCoincidentSides = stk::mesh::impl::get_exposed_coincident_sides(graph, numSidesHex8, 2);
    ASSERT_EQ(0u, exposedCoincidentSides.size());
}

TEST(CoincidentElements, ThreeCoincidentShells)
{
    stk::mesh::Graph graph;
    std::vector<stk::mesh::impl::CoincidentElementDescription> elemDescs = {{2, 0, 1}, {2, 0, 2}, {2, 1, 0}, {2, 1, 2}, {2, 2, 0}, {2, 2, 1}};
    make_graph_for_coincident_elements(graph, 3, elemDescs);
    expect_all_sides_exposed(graph, 2, 0, {1, 2});
    expect_all_sides_exposed(graph, 2, 1, {0, 2});
    expect_all_sides_exposed(graph, 2, 2, {0, 1});
}

void setup_graph_for_hex_next_to_shell(stk::mesh::Graph &graph, const stk::mesh::GraphEdge &graphEdge)
{
    graph.set_num_local_elements(2);
    add_symmetric_edges(graph, graphEdge);
}

TEST(CoincidentElements, Hex8WithShellNeighbor)
{
    stk::mesh::Graph graph;
    stk::mesh::GraphEdge edgeFromHexToShell(0, 1, 1, 0);
    setup_graph_for_hex_next_to_shell(graph, edgeFromHexToShell);

    size_t numSidesHex8 = 6;
    std::vector<stk::mesh::CoincidentElementConnection> exposedCoincidentSides = stk::mesh::impl::get_exposed_coincident_sides(graph, numSidesHex8, edgeFromHexToShell.elem1);
    ASSERT_EQ(0u, exposedCoincidentSides.size());

    size_t numSidesShell = 2;
    exposedCoincidentSides = stk::mesh::impl::get_exposed_coincident_sides(graph, numSidesShell, edgeFromHexToShell.elem2);
    ASSERT_EQ(0u, exposedCoincidentSides.size());
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
