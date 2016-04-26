
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <GameOfLife/NoGhostGameofLife.hpp>  // for NoGhostGameofLife
#include <GameOfLife/PNGProcessor.hpp>  // for PNGProcessor
#include <MeshBuilder/CoordinateSets.hpp>  // for generate_two_dim_elem_id
#include <MeshBuilder/MeshBuilder.hpp>  // for HexMeshBuilder, etc
#include <MeshBuilder/MeshSnake.hpp>    // for HexMeshSnake, QuadMeshSnake, etc
#include <MeshBuilder/MultiImageReader.hpp>  // for MultiImageReader
#include <iostream>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_Wtime, MPI_COMM_WORLD, etc
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityIdVector, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size, etc
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

typedef stk::mesh::Field<int> ScalarIntField;

bool element_ids_are_active(const stk::mesh::BulkData& bulkData, const
                                     stk::mesh::EntityIdVector& elemIds)
{
    bool result = true;
    stk::mesh::Entity elem;
    for (stk::mesh::EntityId elemId : elemIds)
    {
        elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);
        if (!bulkData.is_valid(elem))
            result = false;
    }
    return result;
}

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
        QuadMeshBuilder Mesh1(comm, "1ProcQuadFillArea");
        QuadMeshBuilder Mesh2(comm, "1ProcQuadFillAreaRandomly");
        QuadMeshBuilder Mesh3(comm, "1ProcQuadFillAreaOnProc");
        Mesh1.commit_meta();
        Mesh2.commit_meta();
        Mesh3.commit_meta();

        EXPECT_EQ(0u, Mesh1.num_elems());
        EXPECT_EQ(0u, Mesh2.num_elems());
        EXPECT_EQ(0u, Mesh3.num_elems());

        Mesh1.begin_modification();
        Mesh2.begin_modification();
        Mesh3.begin_modification();
        Mesh1.fill_area(1, 2, 1, 2);
        Mesh2.fill_area_randomly(1,2,1,2);
        Mesh3.fill_area_on_proc(1,2,1,2,0);
        Mesh1.end_modification();
        Mesh2.end_modification();
        Mesh3.end_modification();

        EXPECT_EQ(4u, Mesh1.num_elems());
        EXPECT_EQ(4u, Mesh2.num_elems());
        EXPECT_EQ(4u, Mesh3.num_elems());
    }
}
TEST(MeshBuilder, 4ProcQuadFillArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "4ProcQuadFillArea");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area(1, 5, 1, 5);
        Mesh.end_modification();

        stk::mesh::EntityIdVector elemIds;

        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(2);
            elemIds.push_back(4);
            elemIds.push_back(7);
            elemIds.push_back(11);
            elemIds.push_back(3);
            elemIds.push_back(5);
            elemIds.push_back(8);
            elemIds.push_back(12);
            elemIds.push_back(17);
        }
        else if (1 == procRank)
        {
            elemIds.push_back(6);
            elemIds.push_back(9);
            elemIds.push_back(13);
            elemIds.push_back(18);
            elemIds.push_back(24);
        }
        else if (2 == procRank)
        {
            elemIds.push_back(10);
            elemIds.push_back(14);
            elemIds.push_back(19);
            elemIds.push_back(25);
            elemIds.push_back(32);
        }
        else if (3 == procRank)
        {
            elemIds.push_back(15);
            elemIds.push_back(20);
            elemIds.push_back(26);
            elemIds.push_back(33);
            elemIds.push_back(41);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
    }
}
TEST(MeshBuilder, 4ProcQuadFillAreaRandom)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "4ProcQuadFillAreaRandom");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area_randomly(1, 5, 1, 5);
        Mesh.end_modification();

        stk::mesh::EntityIdVector elemIds;

        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(11);
            elemIds.push_back(12);
            elemIds.push_back(13);
            elemIds.push_back(14);
            elemIds.push_back(15);
            elemIds.push_back(41);
        }
        else if (1 == procRank)
        {
            elemIds.push_back(2);
            elemIds.push_back(3);
            elemIds.push_back(17);
            elemIds.push_back(18);
            elemIds.push_back(19);
            elemIds.push_back(20);
        }
        else if (2 == procRank)
        {
            elemIds.push_back(4);
            elemIds.push_back(5);
            elemIds.push_back(6);
            elemIds.push_back(24);
            elemIds.push_back(25);
            elemIds.push_back(26);
        }
        else if (3 == procRank)
        {
            elemIds.push_back(7);
            elemIds.push_back(8);
            elemIds.push_back(9);
            elemIds.push_back(10);
            elemIds.push_back(32);
            elemIds.push_back(33);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
    }
}
TEST(MeshBuilder, 4ProcQuadFillAreaOnProc)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        QuadMeshBuilder Mesh(comm,"4ProcQuadFillAreaOnProc");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area_on_proc(1, 2, 1, 2, 0);
        Mesh.fill_area_on_proc(1, 2, 3, 4, 1);
        Mesh.fill_area_on_proc(3, 4, 1, 2, 2);
        Mesh.fill_area_on_proc(3, 4, 3, 4, 3);
        Mesh.end_modification();

        EXPECT_EQ(4u, Mesh.num_elems());

        stk::mesh::EntityIdVector elemIds;

        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(2);
            elemIds.push_back(3);
            elemIds.push_back(5);
        }
        else if (1 == procRank)
        {
            elemIds.push_back(6);
            elemIds.push_back(9);
            elemIds.push_back(10);
            elemIds.push_back(14);

        }
        else if (2 == procRank)
        {
            elemIds.push_back(4);
            elemIds.push_back(7);
            elemIds.push_back(8);
            elemIds.push_back(12);

        }

        else if (3 == procRank)
        {
            elemIds.push_back(13);
            elemIds.push_back(18);
            elemIds.push_back(19);
            elemIds.push_back(25);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
    }
}
TEST(MeshBuilder, 4ProcQuadFillAreaLayers)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        QuadMeshBuilder Mesh(comm,"4ProcQuadFillAreaLayers");
        Mesh.commit_meta();

        Mesh.begin_modification();
        Mesh.fill_area_with_layers(1, 4, 1, 6);
        Mesh.end_modification();

        if (0 == procRank || 1 == procRank)
            EXPECT_EQ(8u, Mesh.num_elems());
        else if (2 == procRank || 3 == procRank)
            EXPECT_EQ(4u, Mesh.num_elems());

        stk::mesh::EntityIdVector elemIds;

        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(2);
            elemIds.push_back(4);
            elemIds.push_back(7);
            elemIds.push_back(15);
            elemIds.push_back(20);
            elemIds.push_back(26);
            elemIds.push_back(33);
        }
        else if (1 == procRank)
        {
            elemIds.push_back(3);
            elemIds.push_back(5);
            elemIds.push_back(8);
            elemIds.push_back(12);
            elemIds.push_back(21);
            elemIds.push_back(27);
            elemIds.push_back(34);
            elemIds.push_back(42);
        }
        else if (2 == procRank)
        {
            elemIds.push_back(6);
            elemIds.push_back(9);
            elemIds.push_back(13);
            elemIds.push_back(18);

        }
        else if (3 == procRank)
        {
            elemIds.push_back(10);
            elemIds.push_back(14);
            elemIds.push_back(19);
            elemIds.push_back(25);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
    }
}
TEST(MeshBuilder, 1ProcHexFillArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexMeshBuilder Mesh1(comm, "1ProcHexFillArea");
        HexMeshBuilder Mesh2(comm, "1ProcHexFillAreaOnProc");
        HexMeshBuilder Mesh3(comm, "1ProcHexFillAreaRandom");
        Mesh1.commit_meta();
        Mesh2.commit_meta();
        Mesh3.commit_meta();

        EXPECT_EQ(0u, Mesh1.num_elems());
        EXPECT_EQ(0u, Mesh2.num_elems());
        EXPECT_EQ(0u, Mesh3.num_elems());

        Mesh1.begin_modification();
        Mesh2.begin_modification();
        Mesh3.begin_modification();
        Mesh1.fill_area(1, 2, 1, 2, 1, 2);
        Mesh2.fill_area_randomly(1, 2, 1, 2, 1, 2);
        Mesh3.fill_area_on_proc(1, 2, 1, 2, 1, 2, 0);
        Mesh1.end_modification();
        Mesh2.end_modification();
        Mesh3.end_modification();

        EXPECT_EQ(8u, Mesh1.num_elems());
        EXPECT_EQ(8u, Mesh2.num_elems());
        EXPECT_EQ(8u, Mesh3.num_elems());
    }
}
TEST(MeshBuilder, 4ProcHexFillArea)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "4ProcHexFillAream");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area(1, 2, 1, 2, 1, 5);
        Mesh.end_modification();

        stk::mesh::EntityIdVector elemIds;
        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(2);
            elemIds.push_back(4);
            elemIds.push_back(8);
            elemIds.push_back(3);
            elemIds.push_back(6);
            elemIds.push_back(9);
            elemIds.push_back(16);
        }
        else if (1 == procRank)
        {
            elemIds.push_back(10);
            elemIds.push_back(18);
            elemIds.push_back(19);
            elemIds.push_back(31);
        }
        else if (2 == procRank)
        {
            elemIds.push_back(20);
            elemIds.push_back(33);
            elemIds.push_back(34);
            elemIds.push_back(52);
        }
        else if (3 == procRank)
        {
            elemIds.push_back(35);
            elemIds.push_back(54);
            elemIds.push_back(55);
            elemIds.push_back(80);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
    }
}
TEST(MeshBuilder, 4ProcHexFillAreaRandom)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "4ProcHexFillAreaRandom");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area_randomly(1, 2, 1, 2, 1, 2);
        Mesh.end_modification();

        EXPECT_EQ(2u, Mesh.num_elems());

        stk::mesh::EntityIdVector elemIds;
        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(2);
        }
        else if (1 == procRank)
        {
            elemIds.push_back(4);
            elemIds.push_back(8);
        }
        else if (2 == procRank)
        {
            elemIds.push_back(3);
            elemIds.push_back(6);
        }
        else if (3 == procRank)
        {
            elemIds.push_back(9);
            elemIds.push_back(16);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
    }
}
TEST(MeshBuilder, 4ProcHexFillAreaOnProc)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "4ProcHexFillAreaOnProc");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area_on_proc(1, 2, 1, 1, 1, 1, 0);
        Mesh.fill_area_on_proc(1, 2, 2, 2, 1, 1, 1);
        Mesh.fill_area_on_proc(1, 2, 1, 1, 2, 2, 2);
        Mesh.fill_area_on_proc(1, 2, 2, 2, 2, 2, 3);
        Mesh.end_modification();

        EXPECT_EQ(2u, Mesh.num_elems());

        stk::mesh::EntityIdVector elemIds;
        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(2);
        }
        else if (1 == procRank)
        {
            elemIds.push_back(3);
            elemIds.push_back(6);
        }
        else if (2 == procRank)
        {
            elemIds.push_back(4);
            elemIds.push_back(8);
        }
        else if (3 == procRank)
        {
            elemIds.push_back(9);
            elemIds.push_back(16);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
    }
}
TEST(MeshBuilder, 4ProcHexFillAreaLayers)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexMeshBuilder Mesh(comm, "4ProcHexFillAreaLayers");
        Mesh.commit_meta();

        EXPECT_EQ(0u, Mesh.num_elems());

        Mesh.begin_modification();
        Mesh.fill_area_with_layers(1, 2, 1, 3, 1, 6);
        Mesh.end_modification();

        if (0 == numProcs || 1 == numProcs)
            EXPECT_EQ(12u, Mesh.num_elems());
        else if (2 == numProcs || 3 == numProcs)
            EXPECT_EQ(6u, Mesh.num_elems());

        stk::mesh::EntityIdVector elemIds;
        if (0 == procRank)
        {
            elemIds.push_back(1);
            elemIds.push_back(2);
            elemIds.push_back(3);
            elemIds.push_back(6);
            elemIds.push_back(7);
            elemIds.push_back(13);

            elemIds.push_back(35);
            elemIds.push_back(54);
            elemIds.push_back(55);
            elemIds.push_back(80);
            elemIds.push_back(81);
            elemIds.push_back(113);

        }
        else if (1 == procRank)
        {
            elemIds.push_back(4);
            elemIds.push_back(8);
            elemIds.push_back(9);
            elemIds.push_back(16);
            elemIds.push_back(17);
            elemIds.push_back(28);


            elemIds.push_back(56);
            elemIds.push_back(82);
            elemIds.push_back(83);
            elemIds.push_back(116);
            elemIds.push_back(117);
            elemIds.push_back(158);
        }
        else if (2 == procRank)
        {
            elemIds.push_back(10);
            elemIds.push_back(18);
            elemIds.push_back(19);
            elemIds.push_back(31);
            elemIds.push_back(32);
            elemIds.push_back(49);
        }
        else if (3 == procRank)
        {
            elemIds.push_back(20);
            elemIds.push_back(33);
            elemIds.push_back(34);
            elemIds.push_back(52);
            elemIds.push_back(53);
            elemIds.push_back(77);
        }

        EXPECT_TRUE(element_ids_are_active(Mesh.bulk_data(), elemIds));
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
        Mesh.fill_area_randomly(1, 4, 1, 4);
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
        Mesh.fill_area_randomly(1, 5, 1, 5);
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
        Mesh.fill_area_randomly(1, 4, 1, 4);
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
        Mesh.fill_area_randomly(1, 5, 1, 5);
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
TEST(Mesh3D, DISABLED_FullBodyPerformanceTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procRank = stk::parallel_machine_rank(comm);

    MultiImageReader reader("fullbody", 146);

    {
        HexMeshBuilder meshRandom(comm, "RandomBody");
        meshRandom.commit_meta();
        meshRandom.begin_modification();
        reader.create_randomly_decomposed_mesh(meshRandom);
        double t1 = MPI_Wtime();
        meshRandom.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t1 << std::endl;
        meshRandom.write_mesh();
    }

    {
        HexMeshBuilder meshX(comm, "SingleXBody");
        meshX.commit_meta();
        meshX.begin_modification();
        reader.create_x_layered_decomposed_mesh(meshX);
        double t2 = MPI_Wtime();
        meshX.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t2 << std::endl;
        meshX.write_mesh();
    }

    {
        HexMeshBuilder meshY(comm, "SingleYBody");
        meshY.commit_meta();
        meshY.begin_modification();
        reader.create_y_layered_decomposed_mesh(meshY);
        double t3 = MPI_Wtime();
        meshY.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t3 << std::endl;
        meshY.write_mesh();
    }

    {
        HexMeshBuilder meshZ(comm, "SingleZBody");
        meshZ.commit_meta();
        meshZ.begin_modification();
        reader.create_z_layered_decomposed_mesh(meshZ);
        double t4 = MPI_Wtime();
        meshZ.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t4 << std::endl;
        meshZ.write_mesh();
    }

    {
        HexMeshBuilder MeshX(comm, "BlockXBody");
        MeshX.commit_meta();
        MeshX.begin_modification();
        reader.create_x_blocked_decomposed_mesh(MeshX);
        double t5 = MPI_Wtime();
        MeshX.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t5 << std::endl;
        MeshX.write_mesh();
    }

    {
        HexMeshBuilder MeshY(comm, "BlockYBody");
        MeshY.commit_meta();
        MeshY.begin_modification();
        reader.create_y_blocked_decomposed_mesh(MeshY);
        double t6 = MPI_Wtime();
        MeshY.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t6 << std::endl;
        MeshY.write_mesh();
    }

    {
        HexMeshBuilder MeshZ(comm, "BlockZBody");
        MeshZ.commit_meta();
        MeshZ.begin_modification();
        reader.create_z_blocked_decomposed_mesh(MeshZ);
        double t7 = MPI_Wtime();
        MeshZ.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t7 << std::endl << std::endl;
        MeshZ.write_mesh();
    }
}
TEST(Mesh3D, DISABLED_BroccoliPerformanceTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procRank = stk::parallel_machine_rank(comm);

    MultiImageReader reader("broccoli", 50);

    {
        HexMeshBuilder meshRandom(comm, "RandomBroccoli");
        meshRandom.commit_meta();
        meshRandom.begin_modification();
        reader.create_randomly_decomposed_mesh(meshRandom);
        double t1 = MPI_Wtime();
        meshRandom.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t1 << std::endl;
        meshRandom.write_mesh();
    }

    {
        HexMeshBuilder meshX(comm, "SingleXBroccoli");
        meshX.commit_meta();
        meshX.begin_modification();
        reader.create_x_layered_decomposed_mesh(meshX);
        double t2 = MPI_Wtime();
        meshX.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t2 << std::endl;
        meshX.write_mesh();
    }

    {
        HexMeshBuilder meshY(comm, "SingleYBroccoli");
        meshY.commit_meta();
        meshY.begin_modification();
        reader.create_y_layered_decomposed_mesh(meshY);
        double t3 = MPI_Wtime();
        meshY.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t3 << std::endl;
        meshY.write_mesh();
    }

    {
        HexMeshBuilder meshZ(comm, "SingleZBroccoli");
        meshZ.commit_meta();
        meshZ.begin_modification();
        reader.create_z_layered_decomposed_mesh(meshZ);
        double t4 = MPI_Wtime();
        meshZ.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t4 << std::endl;
        meshZ.write_mesh();
    }

    {
        HexMeshBuilder MeshX(comm, "BlockXBroccoli");
        MeshX.commit_meta();
        MeshX.begin_modification();
        reader.create_x_blocked_decomposed_mesh(MeshX);
        double t5 = MPI_Wtime();
        MeshX.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t5 << std::endl;
        MeshX.write_mesh();
    }

    {
        HexMeshBuilder MeshY(comm, "BlockYBroccoli");
        MeshY.commit_meta();
        MeshY.begin_modification();
        reader.create_y_blocked_decomposed_mesh(MeshY);
        double t6 = MPI_Wtime();
        MeshY.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t6 << std::endl;
        MeshY.write_mesh();
    }

    {
        HexMeshBuilder MeshZ(comm, "BlockZBroccoli");
        MeshZ.commit_meta();
        MeshZ.begin_modification();
        reader.create_z_blocked_decomposed_mesh(MeshZ);
        double t7 = MPI_Wtime();
        MeshZ.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t7 << std::endl << std::endl;
        MeshZ.write_mesh();
    }
}
TEST(Mesh3D, DISABLED_HeadPerformanceTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procRank = stk::parallel_machine_rank(comm);

    MultiImageReader reader("brain3", 80);

    {
        HexMeshBuilder meshRandom(comm, "RandomHead");
        meshRandom.commit_meta();
        meshRandom.begin_modification();
        reader.create_randomly_decomposed_mesh(meshRandom);
        double t1 = MPI_Wtime();
        meshRandom.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t1 << std::endl;
        meshRandom.write_mesh();
    }

    {
        HexMeshBuilder meshX(comm, "SingleXHead");
        meshX.commit_meta();
        meshX.begin_modification();
        reader.create_x_layered_decomposed_mesh(meshX);
        double t2 = MPI_Wtime();
        meshX.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t2 << std::endl;
        meshX.write_mesh();
    }

    {
        HexMeshBuilder meshY(comm, "SingleYHead");
        meshY.commit_meta();
        meshY.begin_modification();
        reader.create_y_layered_decomposed_mesh(meshY);
        double t3 = MPI_Wtime();
        meshY.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t3 << std::endl;
        meshY.write_mesh();
    }

    {
        HexMeshBuilder meshZ(comm, "SingleZHead");
        meshZ.commit_meta();
        meshZ.begin_modification();
        reader.create_z_layered_decomposed_mesh(meshZ);
        double t4 = MPI_Wtime();
        meshZ.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t4 << std::endl;
        meshZ.write_mesh();
    }

    {
        HexMeshBuilder MeshX(comm, "BlockXHead");
        MeshX.commit_meta();
        MeshX.begin_modification();
        reader.create_x_blocked_decomposed_mesh(MeshX);
        double t5 = MPI_Wtime();
        MeshX.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t5 << std::endl;
        MeshX.write_mesh();
    }

    {
        HexMeshBuilder MeshY(comm, "BlockYHead");
        MeshY.commit_meta();
        MeshY.begin_modification();
        reader.create_y_blocked_decomposed_mesh(MeshY);
        double t6 = MPI_Wtime();
        MeshY.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t6 << std::endl;
        MeshY.write_mesh();
    }

    {
        HexMeshBuilder MeshZ(comm, "BlockZHead");
        MeshZ.commit_meta();
        MeshZ.begin_modification();
        reader.create_z_blocked_decomposed_mesh(MeshZ);
        double t7 = MPI_Wtime();
        MeshZ.end_modification();
        if (0 == procRank)
            std::cout << MPI_Wtime() - t7 << std::endl << std::endl;
        MeshZ.write_mesh();
    }
}
TEST(JFF, DISABLED_GoLTuring)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (8 == numProcs)
    {
        QuadMeshBuilder Mesh(comm, "8ProcRainbowGoL");
        stk::mesh::Field<int>* lifeField = nullptr;
        stk::mesh::Field<int>* neighborField = nullptr;
        Mesh.create_life_and_neighbor_fields(lifeField, neighborField);
        Mesh.commit_meta();

        PNGProcessor PNG("turing.png");
        PNG.commit_image_vector_to_pixel_vector();

        Mesh.begin_modification();
        unsigned width = PNG.get_image_width();
        unsigned height = PNG.get_image_height();
        Mesh.fill_area(1, width, 1, height);
        Mesh.end_modification();

        NoGhostGameofLife Game(Mesh.bulk_data(), *lifeField, *neighborField, "RainbowTuring");

        std::vector<std::pair<unsigned, unsigned>> coords;
        PNG.get_coordinates_of_active_pixels(coords);

        stk::mesh::EntityIdVector elemIds;
        for (std::pair<unsigned, unsigned>& pair : coords)
            elemIds.push_back(generate_two_dim_elem_id(pair.first, pair.second));

        Game.activate_these_ids(elemIds);
        Game.run_game_of_life(29);
    }
}
TEST(JFF, DISABLED_Maze)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    QuadMeshBuilder Mesh(comm, "Maze");
    stk::mesh::Field<int>* lifeField = nullptr;
    stk::mesh::Field<int>* neighborField = nullptr;
    Mesh.create_life_and_neighbor_fields(lifeField, neighborField);
    Mesh.commit_meta();

    Mesh.begin_modification();
    Mesh.fill_area(1, 201, 1, 201);
    Mesh.end_modification();

    stk::mesh::EntityIdVector elemIds = {generate_two_dim_elem_id(100, 100)};

    NoGhostGameofLife Game(Mesh.bulk_data(), *lifeField, *neighborField, "Maze");
    Game.activate_these_ids(elemIds);
    Game.run_game_of_life(500);
}
}
