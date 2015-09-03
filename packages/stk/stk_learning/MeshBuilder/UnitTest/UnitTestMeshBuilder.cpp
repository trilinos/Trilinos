#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>

#include <ctime>

#include <GameOfLife/NoGhostGameofLife.hpp>
#include <GameOfLife/PNGProcessor.hpp>
#include <MeshBuilder/MeshBuilder.hpp>
#include <MeshBuilder/MeshSnake.hpp>
#include <MeshBuilder/CoordinateSets.hpp>
#include <MeshBuilder/MultiImageReader.hpp>

typedef stk::mesh::Field<int> ScalarIntField;

namespace
{
TEST(MeshBuilder, 1ProcQuadCreation)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadCreation");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.create_element(1, 1);
        Mesh.end_modification();

        stk::mesh::EntityVector elements;
        stk::mesh::EntityVector nodes;
        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::NODE_RANK, nodes);
        EXPECT_EQ(1u, elements.size());
        EXPECT_EQ(4u, nodes.size());
    }
}
TEST(MeshBuilder, 1ProcQuadNodeCoordinates)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadNodeCoordinates");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.create_element(1, 1);
        Mesh.end_modification();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
        ASSERT_EQ(1u, elements.size());

        stk::mesh::Entity elem = elements[0];

        const stk::mesh::Entity* nodes = Mesh.bulk_data().begin_nodes(elem);
        const stk::mesh::ConnectivityOrdinal* nodeOrds = Mesh.bulk_data().begin_node_ordinals(elem);
        unsigned numNodes = Mesh.bulk_data().num_nodes(elem);

        ASSERT_EQ(4u, numNodes);

        for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
        {
            stk::mesh::Entity node = nodes[nodeIndex];
            stk::mesh::ConnectivityOrdinal nodeOrd = nodeOrds[nodeIndex];

            if (0u == nodeOrd)
            {
                EXPECT_EQ(0u, Mesh.node_x_coord(node));
                EXPECT_EQ(0u, Mesh.node_y_coord(node));
            }
            else if (1u == nodeOrd)
            {
                EXPECT_EQ(1u, Mesh.node_x_coord(node));
                EXPECT_EQ(0u, Mesh.node_y_coord(node));
            }
            else if (2u == nodeOrd)
            {
                EXPECT_EQ(1u, Mesh.node_x_coord(node));
                EXPECT_EQ(1u, Mesh.node_y_coord(node));
            }
            else if (3u == nodeOrd)
            {
                EXPECT_EQ(0u, Mesh.node_x_coord(node));
                EXPECT_EQ(1u, Mesh.node_y_coord(node));
            }
        }
    }
}
TEST(MeshBuilder, 1ProcQuadNodeIds)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadNodeIds");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.create_element(2, 2);
        Mesh.end_modification();

       stk::mesh::EntityVector elements;
       stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
       stk::mesh::Entity elem = elements[0];

       const stk::mesh::Entity* nodes = Mesh.bulk_data().begin_nodes(elem);
       const stk::mesh::ConnectivityOrdinal* nodeOrds = Mesh.bulk_data().begin_node_ordinals(elem);
       unsigned numNodes = Mesh.bulk_data().num_nodes(elem);

       ASSERT_EQ(4u, numNodes);

       for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
       {
           stk::mesh::Entity node = nodes[nodeIndex];
           stk::mesh::ConnectivityOrdinal nodeOrd = nodeOrds[nodeIndex];

           if (0u == nodeOrd)
           {
               EXPECT_EQ(5u, Mesh.bulk_data().identifier(node));
               EXPECT_EQ(1u, Mesh.node_x_coord(node));
               EXPECT_EQ(1u, Mesh.node_y_coord(node));
           }
           else if (1u == nodeOrd)
           {
               EXPECT_EQ(8u, Mesh.bulk_data().identifier(node));
               EXPECT_EQ(2u, Mesh.node_x_coord(node));
               EXPECT_EQ(1u, Mesh.node_y_coord(node));
           }
           else if (2u == nodeOrd)
           {
               EXPECT_EQ(13u, Mesh.bulk_data().identifier(node));
               EXPECT_EQ(2u, Mesh.node_x_coord(node));
               EXPECT_EQ(2u, Mesh.node_y_coord(node));
           }
           else if (3u == nodeOrd)
           {
               EXPECT_EQ(9u, Mesh.bulk_data().identifier(node));
               EXPECT_EQ(1u, Mesh.node_x_coord(node));
               EXPECT_EQ(2u, Mesh.node_y_coord(node));
           }
       }
    }
}
TEST(MeshBuilder, 4ProcShareNodes)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == procSize)
    {
       QuadMeshBuilder Mesh(comm, "4ProcShareNodes");
       Mesh.commit_meta();

       Mesh.begin_modification();
       Mesh.create_element(1, 1, 0);
       Mesh.create_element(1, 2, 1);
       Mesh.create_element(2, 1, 2);
       Mesh.create_element(2, 2, 3);
       Mesh.end_modification();

       EXPECT_EQ(1u, Mesh.num_elems());

       stk::mesh::Entity node5 = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 5);
       std::vector<int> node5procs;
       Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(node5), node5procs);
       EXPECT_EQ(3u, node5procs.size());

       if (1 == procRank || 0 == procRank)
       {
           stk::mesh::Entity node3 = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 3);
           std::vector<int> node3procs;
           Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(node3), node3procs);
           EXPECT_EQ(1u, node3procs.size());
       }
       else if (3 == procRank || 2 == procRank)
       {
           stk::mesh::Entity node8 = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 8);
           std::vector<int> node8procs;
           Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(node8), node8procs);
           EXPECT_EQ(1u, node8procs.size());
       }
    }
}
//TEST(MeshBuilder, 4ProcOwnNodes)
//{
//    stk::ParallelMachine comm = MPI_COMM_WORLD;
//    int procSize = stk::parallel_machine_size(comm);
//    int procRank = stk::parallel_machine_rank(comm);
//    if (4 == procSize)
//    {
//       QuadMeshBuilder Mesh(comm, "4ProcOwnNodes");
//       Mesh.commit_meta();
//
//       Mesh.begin_modification();
//       Mesh.create_element(1, 1);
//       Mesh.create_element(2, 1);
//       Mesh.create_element(1, 2);
//       Mesh.create_element(2, 2);
//       Mesh.end_modification();
//
//
//
//    }
//}
TEST(MeshBuilder, 1ProcHexCreation)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        HexMeshBuilder Mesh(comm, "1ProcHexCreation");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.create_element(1, 1, 1);
        Mesh.end_modification();

        stk::mesh::EntityVector elements;
        stk::mesh::EntityVector nodes;
        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::NODE_RANK, nodes);
        EXPECT_EQ(1u, elements.size());
        EXPECT_EQ(8u, nodes.size());
    }
}
TEST(MeshBuilder, 1ProcHexNodeCoordinates)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        HexMeshBuilder Mesh(comm, "1ProcHexNodeCoordinates");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.create_element(1, 1, 1);
        Mesh.end_modification();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
        ASSERT_EQ(1u, elements.size());

        stk::mesh::Entity elem = elements[0];

        const stk::mesh::Entity* nodes = Mesh.bulk_data().begin_nodes(elem);
        const stk::mesh::ConnectivityOrdinal* nodeOrds = Mesh.bulk_data().begin_node_ordinals(elem);
        unsigned numNodes = Mesh.bulk_data().num_nodes(elem);

        ASSERT_EQ(8u, numNodes);

        for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
        {
            stk::mesh::Entity node = nodes[nodeIndex];
            stk::mesh::ConnectivityOrdinal nodeOrd = nodeOrds[nodeIndex];

            if (0u == nodeOrd)
            {
                EXPECT_EQ(0u, Mesh.node_x_coord(node));
                EXPECT_EQ(0u, Mesh.node_y_coord(node));
                EXPECT_EQ(0u, Mesh.node_z_coord(node));
            }
            else if (1u == nodeOrd)
            {
                EXPECT_EQ(1u, Mesh.node_x_coord(node));
                EXPECT_EQ(0u, Mesh.node_y_coord(node));
                EXPECT_EQ(0u, Mesh.node_z_coord(node));
            }
            else if (2u == nodeOrd)
            {
                EXPECT_EQ(1u, Mesh.node_x_coord(node));
                EXPECT_EQ(1u, Mesh.node_y_coord(node));
                EXPECT_EQ(0u, Mesh.node_z_coord(node));
            }
            else if (3u == nodeOrd)
            {
                EXPECT_EQ(0u, Mesh.node_x_coord(node));
                EXPECT_EQ(1u, Mesh.node_y_coord(node));
                EXPECT_EQ(0u, Mesh.node_z_coord(node));
            }
            else if (4u == nodeOrd)
            {
                EXPECT_EQ(0u, Mesh.node_x_coord(node));
                EXPECT_EQ(0u, Mesh.node_y_coord(node));
                EXPECT_EQ(1u, Mesh.node_z_coord(node));
            }
            else if (5u == nodeOrd)
            {
                EXPECT_EQ(1u, Mesh.node_x_coord(node));
                EXPECT_EQ(0u, Mesh.node_y_coord(node));
                EXPECT_EQ(1u, Mesh.node_z_coord(node));
            }
            else if (6u == nodeOrd)
            {
                EXPECT_EQ(1u, Mesh.node_x_coord(node));
                EXPECT_EQ(1u, Mesh.node_y_coord(node));
                EXPECT_EQ(1u, Mesh.node_z_coord(node));
            }
            else if (7u == nodeOrd)
            {
                EXPECT_EQ(0u, Mesh.node_x_coord(node));
                EXPECT_EQ(1u, Mesh.node_y_coord(node));
                EXPECT_EQ(1u, Mesh.node_z_coord(node));
            }
        }
    }
}
TEST(MeshBuilder, 1ProcHexNodeIds)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        HexMeshBuilder Mesh(comm, "1ProcHexNodeIds");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.create_element(1,1,1);
        Mesh.end_modification();

       stk::mesh::EntityVector elements;
       stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
       stk::mesh::Entity elem = elements[0];

       const stk::mesh::Entity* nodes = Mesh.bulk_data().begin_nodes(elem);
       const stk::mesh::ConnectivityOrdinal* nodeOrds = Mesh.bulk_data().begin_node_ordinals(elem);
       unsigned numNodes = Mesh.bulk_data().num_nodes(elem);

       ASSERT_EQ(8u, numNodes);

       for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
       {
           stk::mesh::Entity node = nodes[nodeIndex];
           stk::mesh::ConnectivityOrdinal nodeOrd = nodeOrds[nodeIndex];

           if (0u == nodeOrd)
           {
               EXPECT_EQ(1u, Mesh.bulk_data().identifier(node));
           }
           else if (1u == nodeOrd)
           {
               EXPECT_EQ(2u, Mesh.bulk_data().identifier(node));
           }
           else if (2u == nodeOrd)
           {
               EXPECT_EQ(6u, Mesh.bulk_data().identifier(node));
           }
           else if (3u == nodeOrd)
           {
               EXPECT_EQ(3u, Mesh.bulk_data().identifier(node));
           }
           else if (4u == nodeOrd)
           {
               EXPECT_EQ(4u, Mesh.bulk_data().identifier(node));
           }
           else if (5u == nodeOrd)
           {
               EXPECT_EQ(8u, Mesh.bulk_data().identifier(node));
           }
           else if (6u == nodeOrd)
           {
               EXPECT_EQ(16u, Mesh.bulk_data().identifier(node));
           }
           else if (7u == nodeOrd)
           {
               EXPECT_EQ(9u, Mesh.bulk_data().identifier(node));
           }
       }
    }
}
TEST(MeshBuilder, 4ProcHexShareNodes)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == procSize)
    {
       HexMeshBuilder Mesh(comm, "4ProcHexShareNodes");
       Mesh.commit_meta();

       Mesh.begin_modification();
       Mesh.create_element(1,1,1,0);
       Mesh.create_element(1,1,2,0);
       Mesh.create_element(2,1,1,1);
       Mesh.create_element(2,1,2,1);
       Mesh.create_element(1,2,1,2);
       Mesh.create_element(1,2,2,2);
       Mesh.create_element(2,2,1,3);
       Mesh.create_element(2,2,2,3);
       Mesh.end_modification();

       EXPECT_EQ(2u, Mesh.num_elems());

       stk::mesh::Entity centerNode = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 16);
       std::vector<int> centerNodeProcs;
       Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(centerNode), centerNodeProcs);
       EXPECT_EQ(3u, centerNodeProcs.size());

       if (0 == procRank || 1 == procRank)
       {
           stk::mesh::Entity node8 = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 8);
           std::vector<int> node8Procs;
           Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(node8), node8Procs);
           EXPECT_EQ(1u, node8Procs.size());
       }
       if (0 == procRank || 2 == procRank)
       {
           stk::mesh::Entity node9 = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 9);
           std::vector<int> node9Procs;
           Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(node9), node9Procs);
           EXPECT_EQ(1u, node9Procs.size());
       }
       if (1 == procRank || 3 == procRank)
       {
           stk::mesh::Entity node27 = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 27);
           std::vector<int> node27Procs;
           Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(node27), node27Procs);
           EXPECT_EQ(1u, node27Procs.size());
       }
       if (2 == procRank || 3 == procRank)
       {
           stk::mesh::Entity node28 = Mesh.bulk_data().get_entity(stk::topology::NODE_RANK, 28);
           std::vector<int> node28Procs;
           Mesh.bulk_data().comm_shared_procs(Mesh.bulk_data().entity_key(node28), node28Procs);
           EXPECT_EQ(1u, node28Procs.size());
       }
    }
}
TEST(MeshBuilder, 1ProcQuadFillArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadFillArea");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area(1, 2, 1, 2);
        Mesh.end_modification();

        EXPECT_EQ(4u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 4ProcQuadFillArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (4 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "4ProcQuadFillArea");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area(1, 4, 1, 4);
        Mesh.end_modification();

        EXPECT_EQ(4u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 1ProcHexFillArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "1ProcHexFillArea");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area(1, 2, 1, 2, 1, 2);
        Mesh.end_modification();

        EXPECT_EQ(8u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 4ProcHexFillArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (4 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "4ProcHexFillArea");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area(1, 4, 1, 4, 1, 4);
        Mesh.end_modification();

        EXPECT_EQ(16u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 1ProcQuadRemove)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadRemove");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 2, 1, 2);
        Mesh.end_modification();

        EXPECT_EQ(4u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_element(1, 1);
        Mesh.begin_modification();

        EXPECT_EQ(3u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 4ProcQuadRemove)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "4ProcQuadRemove");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 4, 1, 4);
        Mesh.end_modification();

        EXPECT_EQ(4u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_element(1, 1);
        Mesh.end_modification();

        if (0 == procRank)
            EXPECT_EQ(3u, Mesh.num_elems());
        else
            EXPECT_EQ(4u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_element(1, 2);
        Mesh.remove_element(1, 3);
        Mesh.remove_element(1, 4);
        Mesh.end_modification();

        EXPECT_EQ(3u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 1ProcHexRemove)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "1ProcHexRemove");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 2, 1, 2, 1, 2);
        Mesh.end_modification();

        EXPECT_EQ(8u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_element(1, 1, 1);
        Mesh.begin_modification();

        EXPECT_EQ(7u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 4ProcHexRemove)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "4ProcHexRemove");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 4, 1, 4, 1, 4);
        Mesh.end_modification();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
        EXPECT_EQ(16u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_element(1, 1, 1);
        Mesh.end_modification();

        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);

        if (0 == procRank)
            EXPECT_EQ(15u, Mesh.num_elems());
        else
            EXPECT_EQ(16u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_element(1, 1, 2);
        Mesh.remove_element(1, 1, 3);
        Mesh.remove_element(1, 1, 4);
        Mesh.end_modification();

        stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
        EXPECT_EQ(15u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 1ProcQuadRemoveArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadRemoveArea");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 5, 1, 5);
        Mesh.end_modification();

        EXPECT_EQ(25u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_area(1, 3, 1, 3);
        Mesh.end_modification();

        EXPECT_EQ(16u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 4ProcQuadRemoveArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "4ProcQuadRemoveArea");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 4, 1, 4);
        Mesh.end_modification();

        EXPECT_EQ(4u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_area(1, 3, 1, 3);
        Mesh.end_modification();

        if (3 == procRank)
            EXPECT_EQ(4u, Mesh.num_elems());
        else
            EXPECT_EQ(1u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 1ProcHexRemoveArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "1ProcHexRemoveArea");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 3, 1, 3, 1, 3);
        Mesh.end_modification();

        EXPECT_EQ(27u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_area(1, 2, 1, 2, 1, 2);
        Mesh.end_modification();

        EXPECT_EQ(19u, Mesh.num_elems());
    }
}
TEST(MeshBuilder, 4ProcHexRemoveArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (4 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "4ProcHexRemoveArea");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 4, 1, 4, 1, 4);
        Mesh.end_modification();

        EXPECT_EQ(16u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.remove_area(1, 4, 1, 1, 1, 4);
        Mesh.end_modification();

        EXPECT_EQ(12u, Mesh.num_elems());
    }
}
TEST(MeshSnake, 1ProcQuadSnakeConstruction)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadSnakeConstruction");
        Mesh.commit_meta();

        QuadMeshSnake Snake(Mesh);
        EXPECT_EQ(1u, Snake.x_lower_bound());
        EXPECT_EQ(1u, Snake.x_upper_bound());
        EXPECT_EQ(1u, Snake.y_lower_bound());
        EXPECT_EQ(1u, Snake.y_upper_bound());
        EXPECT_EQ(1u, Snake.x_pos());
        EXPECT_EQ(1u, Snake.y_pos());
    }
}
TEST(MeshSnake, 1ProcQuadSnakeBounds)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadSnakeBounds");
        Mesh.commit_meta();

        QuadMeshSnake Snake(Mesh);
        Snake.set_x_bounds(2, 10);
        Snake.set_y_bounds(1, 5);
        Snake.set_x_pos(3);
        Snake.set_y_pos(2);

        EXPECT_EQ(2u, Snake.x_lower_bound());
        EXPECT_EQ(10u, Snake.x_upper_bound());
        EXPECT_EQ(1u, Snake.y_lower_bound());
        EXPECT_EQ(5u, Snake.y_upper_bound());
        EXPECT_EQ(3u, Snake.x_pos());
        EXPECT_EQ(2u, Snake.y_pos());

        Snake.set_x_pos(1);
        Snake.set_y_pos(10);

        EXPECT_EQ(3u, Snake.x_pos());
        EXPECT_EQ(2u, Snake.y_pos());

        Snake.set_pos(6, 5);

        EXPECT_EQ(6u, Snake.x_pos());
        EXPECT_EQ(5u, Snake.y_pos());

    }
}
TEST(MeshSnake, 1ProcQuadSnakeBegin)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "1ProcQuadSnakeBegin");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 5, 1, 5);
        Mesh.end_modification();

        QuadMeshSnake Snake(Mesh);
        Snake.set_x_bounds(1, 5);
        Snake.set_y_bounds(1, 5);
        Snake.set_x_pos(1);
        Snake.set_y_pos(1);

        EXPECT_EQ(INVALID_DIR, Snake.dir());

        Snake.begin_snake();

        EXPECT_NE(INVALID_DIR, Snake.dir());
    }
}
TEST(MeshSnake, 1ProcHexSnakeConstruction)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "1ProcQuadSnakeConstruction");
        Mesh.commit_meta();

        HexMeshSnake Snake(Mesh);
        EXPECT_EQ(1u, Snake.x_lower_bound());
        EXPECT_EQ(1u, Snake.x_upper_bound());
        EXPECT_EQ(1u, Snake.y_lower_bound());
        EXPECT_EQ(1u, Snake.y_upper_bound());
        EXPECT_EQ(1u, Snake.z_lower_bound());
        EXPECT_EQ(1u, Snake.z_upper_bound());
        EXPECT_EQ(1u, Snake.x_pos());
        EXPECT_EQ(1u, Snake.y_pos());
        EXPECT_EQ(1u, Snake.z_pos());
    }
}
TEST(MeshSnake, 1ProcHexSnakeBounds)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "1ProcHexSnakeBounds");
        Mesh.commit_meta();

        HexMeshSnake Snake(Mesh);
        Snake.set_x_bounds(2, 10);
        Snake.set_y_bounds(1, 5);
        Snake.set_z_bounds(4, 12);
        Snake.set_x_pos(3);
        Snake.set_y_pos(2);
        Snake.set_z_pos(5);

        EXPECT_EQ(2u, Snake.x_lower_bound());
        EXPECT_EQ(10u, Snake.x_upper_bound());

        EXPECT_EQ(1u, Snake.y_lower_bound());
        EXPECT_EQ(5u, Snake.y_upper_bound());

        EXPECT_EQ(4u, Snake.z_lower_bound());
        EXPECT_EQ(12u, Snake.z_upper_bound());

        EXPECT_EQ(3u, Snake.x_pos());
        EXPECT_EQ(2u, Snake.y_pos());
        EXPECT_EQ(5u, Snake.z_pos());

        Snake.set_x_pos(1);
        Snake.set_y_pos(10);
        Snake.set_z_pos(1);

        EXPECT_EQ(3u, Snake.x_pos());
        EXPECT_EQ(2u, Snake.y_pos());
        EXPECT_EQ(5u, Snake.z_pos());

        Snake.set_pos(6, 5, 6);

        EXPECT_EQ(6u, Snake.x_pos());
        EXPECT_EQ(5u, Snake.y_pos());
        EXPECT_EQ(6u, Snake.z_pos());
    }
}
TEST(MeshSnake, 1ProcHexSnakeBegin)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "1ProcQuadSnakeBegin");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 5, 1, 5, 1, 5);
        Mesh.end_modification();

        HexMeshSnake Snake(Mesh);
        Snake.set_x_bounds(1, 5);
        Snake.set_y_bounds(1, 5);
        Snake.set_z_bounds(1, 5);
        Snake.set_x_pos(2);
        Snake.set_y_pos(2);
        Snake.set_z_pos(2);

        EXPECT_EQ(INVALID_DIR, Snake.dir());

        Snake.begin_snake();

        EXPECT_NE(INVALID_DIR, Snake.dir());
    }
}
TEST(MeshSnake, 8ProcHexSnakeCrawl)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (8 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "8ProcHexSnakeCrawl");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 25, 1, 25, 1, 25);
        Mesh.end_modification();

        HexMeshSnake Snake(Mesh);
        Snake.set_x_bounds(1, 25);
        Snake.set_y_bounds(1, 25);
        Snake.set_z_bounds(1, 25);

        Snake.begin_snake();

        Snake.crawl(10000);
    }
}
//TEST(Mesh3D, MultiImage)
//{
//    MultiImageReader Image("melon", 90);
//}

TEST(Mesh3D, DISABLED_Body)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (8 == numProcs)
    {
        unsigned numFiles = 90;

        HexMeshBuilder Mesh(comm, "Melon");

        unsigned procCounter = 0;
        std::vector<std::pair<unsigned, unsigned>> sliceCoordinates;

        Mesh.begin_modification();
        for (unsigned index = 0; index < numFiles; index++)
        {
            sliceCoordinates.clear();

            std::string fileName = "melon"+std::to_string(index)+".png";
            ColoredPNGProcessor image(fileName);
            image.commit_image_vector_to_pixel_vector_with_greyscale();

            image.get_coordinates_of_inactive_pixels(sliceCoordinates);

            for (std::pair<unsigned, unsigned> xyCoord : sliceCoordinates)
                Mesh.create_element(xyCoord.first, xyCoord.second, index+1, procCounter%numProcs);

        Mesh.write_mesh();
            std::cout << index << std::endl;
            procCounter++;
        }

        std::cout << "done" << std::endl;

        Mesh.end_modification();
    }
}
TEST(MeshSnake, DISABLED_8ProcCrawl)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (8 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "8ProcQuadSnakeCrawl");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1, 100, 1, 100);
        Mesh.end_modification();

        QuadMeshSnake Snake(Mesh);
        Snake.set_x_bounds(1, 100);
        Snake.set_y_bounds(1, 100);
        Snake.set_x_pos(1);
        Snake.set_y_pos(1);
        Snake.begin_snake();

        Snake.crawl(10000);
    }
}
TEST(JFF, DISABLED_GoL)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (8 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "8ProcRainbowGoL");

        stk::mesh::Field<int>* lifeField = &Mesh.meta_data().declare_field<stk::mesh::Field<int>>(
                stk::topology::ELEM_RANK, "lifeField");
        stk::mesh::Field<int>* neighborField = &Mesh.meta_data().declare_field<stk::mesh::Field<int>>(
                stk::topology::ELEM_RANK, "neighborField");
        int val = 0;
        stk::mesh::put_field(*lifeField, Mesh.meta_data().universal_part(), &val);
        stk::mesh::put_field(*neighborField, Mesh.meta_data().universal_part(), &val);

        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area(1,100, 1, 100);
        Mesh.end_modification();

        NoGhostGameofLife Game(&Mesh.bulk_data(), lifeField, neighborField, "Rainbow");

        stk::mesh::EntityIdVector elemIds;

        for (unsigned x = 1; x <= 100; x++)
            for (unsigned y = 1; y <= 100; y++)
                if (!((100*x + y)%5) || !((100*x + y)%7))
                    elemIds.push_back(generate_two_dim_elem_id(x, y));

        Game.activate_these_ids(elemIds);
        Game.run_game_of_life(100);
    }
}
TEST(JFF, DISABLED_PNGPrint)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    QuadMeshBuilder Mesh(comm, "RainbowOceanPrint");
    Mesh.commit_meta();

    ColoredPNGProcessor PNG("ocean.png");
    PNG.commit_image_vector_to_pixel_vector();
    std::vector<std::pair<unsigned, unsigned>> activePoints;
    PNG.get_coordinates_of_active_pixels(activePoints);

    std::random_shuffle(activePoints.begin(), activePoints.end());

    Mesh.begin_modification();
    for (unsigned index = 0, size = activePoints.size(); index < size; index++)
    {
        Mesh.create_element(activePoints[index].first, activePoints[index].second);
        if (!(index%500))
            Mesh.write_mesh();

        std::cout << index << std::endl;
    }
    Mesh.end_modification();
    Mesh.write_mesh();
}
TEST(JFF, DISABLED_PNGErase)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    QuadMeshBuilder Mesh(comm, "RainbowSandiaErase");
    Mesh.commit_meta();

    ColoredPNGProcessor PNG("index.png");
    PNG.commit_image_vector_to_pixel_vector();
    std::vector<std::pair<unsigned, unsigned>> inactivePoints;
    PNG.get_coordinates_of_inactive_pixels(inactivePoints);

    std::random_shuffle(inactivePoints.begin(), inactivePoints.end());

    Mesh.begin_modification();
    Mesh.fill_area(1, PNG.get_image_width(), 1, PNG.get_image_height());
    Mesh.end_modification();

    Mesh.begin_modification();
    for (unsigned index = 0, size = inactivePoints.size(); index < size; index++)
    {
        Mesh.remove_element(inactivePoints[index].first, inactivePoints[index].second);
        if (!(index%100))
            Mesh.write_mesh();

        std::cout << index << std::endl;
    }
    Mesh.end_modification();
    Mesh.write_mesh();
}
}
