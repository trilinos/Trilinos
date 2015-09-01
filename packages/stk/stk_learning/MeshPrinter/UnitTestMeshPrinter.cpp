#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_util/stk_util/parallel/CommSparse.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>

#include <ctime>
#include "MeshPrinter.h"

namespace
{
TEST(MeshPrinter, 1ProcQuadCreation)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        QuadMeshPrinter Mesh(comm, "mesh");
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

TEST(MeshPrinter, 1ProcQuadNodeCoordinates)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        QuadMeshPrinter Mesh(comm, "mesh");
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
TEST(MeshPrinter, 1ProcQuadNodeIds)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        QuadMeshPrinter Mesh(comm, "mesh");
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
TEST(MeshPrinter, 4ProcShareNodes)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == procSize)
    {
       QuadMeshPrinter Mesh(comm, "4ProcMesh");
       Mesh.commit_meta();

       Mesh.begin_modification();
       Mesh.create_element(1, 1, 0);
       Mesh.create_element(1, 2, 1);
       Mesh.create_element(2, 1, 2);
       Mesh.create_element(2, 2, 3);
       Mesh.end_modification();

       Mesh.write_mesh();

       stk::mesh::EntityVector elements;
       stk::mesh::get_entities(Mesh.bulk_data(), stk::topology::ELEM_RANK, elements);
       EXPECT_EQ(1u, elements.size());

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

TEST(MeshPrinter, 8ProcDesign)
{

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);
    if (8 == procSize)
    {
        QuadMeshPrinter Mesh(comm, "8ProcSnake");
        Mesh.commit_meta();

        srand(time(NULL));
        Mesh.begin_modification();
        unsigned x = 500;
        unsigned y = 500;
        for (unsigned step = 1; step <= 1000; step++)
        {
            int procNum = rand()%procSize;
            Mesh.create_element(x, y, procNum);
            Mesh.write_mesh();
            int direction = rand()%3;
            switch (direction)
            {
                case 0:
                    x--;
                    break;
                case 1:
                    x++;
                    break;
                case 2:
                    y--;
                    break;
            }
        }
        Mesh.end_modification();
    }

}
TEST(MeshPrinter, 1ProcHexCreation)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        HexMeshPrinter Mesh(comm, "mesh");
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
TEST(MeshPrinter, 1ProcHexNodeCoordinates)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        HexMeshPrinter Mesh(comm, "mesh");
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
TEST(MeshPrinter, 1ProcHexNodeIds)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procSize = stk::parallel_machine_size(comm);

    if (1 == procSize)
    {
        HexMeshPrinter Mesh(comm, "mesh");
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

}
