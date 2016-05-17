
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <GameOfLife/GameofLife.hpp>    // for PartGameofLife, etc
#include <GameOfLife/GameofLifeMesh.hpp>  // for QuadGameofLifeMesh, etc
#include <GameOfLife/NoGhostGameofLife.hpp>  // for NoGhostGameofLife
#include <GameOfLife/PNGProcessor.hpp>  // for PNGProcessor, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityIdVector, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size, etc
#include "../../../stk_unit_test_utils/getOption.h"





/*
 * There are 1 proc and 4 proc tests, mostly clones of each other to test parallel
 * consistency.
 */
namespace
{
TEST(GameofLifeClass, 1ProcTestTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string meshName = "1ProcTestTest";
        QuadGameofLifeMesh Mesh(comm, 4, 4);

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16};
        PartGameofLife PartGame(Mesh, meshName);
        PartGame.activate_these_ids(elemIds);
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));

        stk::mesh::EntityIdVector notElemIds1 = {1, 2, 3, 4, 5, 6, 8, 9, 12, 13, 14, 15, 16};
        stk::mesh::EntityIdVector notElemIds2 = {1, 2, 3, 4, 5, 8, 9, 10, 12, 13, 14, 15, 16};
        stk::mesh::EntityIdVector elemIds1 = {1, 2, 3, 4};
        stk::mesh::EntityIdVector elemIds2 = {13, 14, 15, 16};

        EXPECT_EQ(16u, PartGame.get_num_elems_on_proc());
        EXPECT_EQ(12u, PartGame.get_num_active_elements());
        EXPECT_FALSE(PartGame.are_these_ids_active(notElemIds1));
        EXPECT_FALSE(PartGame.are_these_ids_active(notElemIds2));
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds1));
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds2));
    }
}
TEST(GameofLifeClass, 4ProcTestTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string meshName = "4ProcTestTest";
        QuadGameofLifeMesh Mesh(comm, 4, 4);

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16};
        PartGameofLife PartGame(Mesh, meshName);
        PartGame.activate_these_ids(elemIds);
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));

        stk::mesh::EntityIdVector notElemIds1 = {1, 2, 3, 4, 5, 6, 8, 9, 12, 13, 14, 15, 16};
        stk::mesh::EntityIdVector notElemIds2 = {1, 2, 3, 4, 5, 8, 9, 10, 12, 13, 14, 15, 16};
        stk::mesh::EntityIdVector elemIds1 = {1, 2, 3, 4};
        stk::mesh::EntityIdVector elemIds2 = {13, 14, 15, 16};
        EXPECT_EQ(4u, PartGame.get_num_elems_on_proc());
        if (1 == procRank || 2 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
        }
        else if (0 == procRank || 3 == procRank)
        {
            EXPECT_EQ(4u, PartGame.get_num_active_elements());
        }
        if (1 == procRank)
        {
            EXPECT_FALSE(PartGame.are_these_ids_active(notElemIds1));
        }
        if (2 == procRank)
        {
            EXPECT_FALSE(PartGame.are_these_ids_active(notElemIds2));
        }
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds1));
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds2));
    }
}
TEST(TriangleGameofLifeClass, 1ProcRandomTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string partMeshName = "1ProcPartRandomTest";
        std::string fieldMeshName = "1ProcFieldRandomTest";

        TriGameofLifeMesh PartMesh(comm, 4, 4);
        TriGameofLifeMesh FieldMesh(comm, 4, 4);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, partMeshName);

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4, 5, 6, 7, 8,
                                             25, 26, 27, 28, 29, 30, 31, 32};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);
        EXPECT_EQ(16u, PartGame.get_num_active_elements());
        EXPECT_EQ(16u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        EXPECT_EQ(4u, PartGame.get_num_active_elements());
        EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        stk::mesh::EntityIdVector exElemIds = {12, 14, 19, 21};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        EXPECT_EQ(4u, PartGame.get_num_active_elements());
        EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
    }
}
TEST(TriangleGameofLifeClass, 4ProcRandomTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string partMeshName = "4ProcPartRandomTest";
        std::string fieldMeshName = "4ProcFieldRandomTest";

        TriGameofLifeMesh PartMesh(comm, 4, 4);
        TriGameofLifeMesh FieldMesh(comm, 4, 4);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, fieldMeshName);

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4, 5, 6, 7, 8,
                                             25, 26, 27, 28, 29, 30, 31, 32};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);
        if (0 == procRank || 3 == procRank)
        {
            EXPECT_EQ(8u, PartGame.get_num_active_elements());
            EXPECT_EQ(8u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);
        if (1 == procRank || 2 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
            EXPECT_EQ(2u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());

        }
        stk::mesh::EntityIdVector exElemIds = {12, 14, 19, 21};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);
        if (1 == procRank || 2 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
            EXPECT_EQ(2u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());

        }
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
    }
}
TEST(TriangleGameofLifeClass, 1ProcInfiniteTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string partMeshName = "1ProcPartInfiniteTest";
        std::string fieldMeshName = "1ProcFieldInfiniteTest";

        TriGameofLifeMesh PartMesh(comm, 10, 5);
        TriGameofLifeMesh FieldMesh(comm, 10, 5);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, fieldMeshName);

        EXPECT_EQ(100u, PartGame.get_num_elems_on_proc());
        EXPECT_EQ(100u, FieldGame.get_num_elems_on_proc());

        stk::mesh::EntityIdVector elemIds = {32, 48, 67, 68};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);
        EXPECT_EQ(4u, PartGame.get_num_active_elements());
        EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds));

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        stk::mesh::EntityIdVector exElemIds1 = {8, 25, 29, 43, 48, 63, 67, 83};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds1));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds1));
        EXPECT_EQ(8u, PartGame.get_num_active_elements());
        EXPECT_EQ(8u, FieldGame.get_num_active_elements());


        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        stk::mesh::EntityIdVector exElemIds2 = {42, 43, 44, 45, 47, 48, 66, 67, 68, 81, 82,
                                                84, 85, 86, 87};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds2));
    }
}
TEST(TriangleGameofLifeClass, 4ProcInfiniteTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string partMeshName = "4ProcPartInfiniteTest";
        std::string fieldMeshName = "4ProcFieldInfiniteTest";

        TriGameofLifeMesh PartMesh(comm, 10, 5);
        TriGameofLifeMesh FieldMesh(comm, 10, 5);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, fieldMeshName);

        if (3 == procRank)
        {
            EXPECT_EQ(40u, PartGame.get_num_elems_on_proc());
            EXPECT_EQ(40u, FieldGame.get_num_elems_on_proc());
        }
        else
        {
            EXPECT_EQ(20u, PartGame.get_num_elems_on_proc());
            EXPECT_EQ(20u, FieldGame.get_num_elems_on_proc());
        }

        stk::mesh::EntityIdVector elemIds = {32, 48, 67, 68};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);
        if (3 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
            EXPECT_EQ(2u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }
        else if (0 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        stk::mesh::EntityIdVector exElemIds1 = {8, 25, 29, 43, 48, 63, 67, 83};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds1));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds1));
        if (3 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
            EXPECT_EQ(2u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
            EXPECT_EQ(2u, FieldGame.get_num_active_elements());
        }
        else if (0 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        stk::mesh::EntityIdVector exElemIds2 = {42, 43, 44, 45, 47, 48, 66, 67, 68,
                                                81, 82, 84, 85, 86, 87};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds2));
        if (3 == procRank)
        {
            EXPECT_EQ(9u, PartGame.get_num_active_elements());
            EXPECT_EQ(9u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(6u, PartGame.get_num_active_elements());
            EXPECT_EQ(6u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
    }
}
TEST(QuadGameofLifeClass, 1ProcGliderTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string meshName1 = "1ProcPartGliderTest";
        std::string meshName2 = "1ProcFieldGliderTest";

        QuadGameofLifeMesh PartMesh(comm, 10, 10);
        QuadGameofLifeMesh FieldMesh(comm, 10, 10);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName1);

        stk::mesh::EntityIdVector elemIds = {71, 72, 73, 83, 92};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        EXPECT_EQ(5u, PartGame.get_num_active_elements());
        EXPECT_EQ(5u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(10);
        FieldGame.run_game_of_life(10);

        EXPECT_EQ(5u, PartGame.get_num_active_elements());
        EXPECT_EQ(5u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(10);
        FieldGame.run_game_of_life(10);

        EXPECT_EQ(5u, PartGame.get_num_active_elements());
        EXPECT_EQ(5u, FieldGame.get_num_active_elements());
    }
}
TEST(QuadGameofLifeClass, 4ProcGliderTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string meshName1 = "4ProcPartGliderTest";
        std::string meshName2 = "4ProcFieldGliderTest";

        QuadGameofLifeMesh PartMesh(comm, 10, 10);
        QuadGameofLifeMesh FieldMesh(comm, 10, 10);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        if (3 == procRank)
        {
            EXPECT_EQ(40u, PartGame.get_num_elems_on_proc());
            EXPECT_EQ(40u, FieldGame.get_num_elems_on_proc());
        }
        else
        {
            EXPECT_EQ(20u, PartGame.get_num_elems_on_proc());
            EXPECT_EQ(20u, FieldGame.get_num_elems_on_proc());
        }

        stk::mesh::EntityIdVector elemIds = {71, 72, 73, 83, 92};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);
        if (3 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
        PartGame.run_game_of_life(10);
        FieldGame.run_game_of_life(10);
        if (3 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(4u, PartGame.get_num_active_elements());
            EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
        PartGame.run_game_of_life(10);
        FieldGame.run_game_of_life(10);
        if (2 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(4u, PartGame.get_num_active_elements());
            EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
    }
}
TEST(QuadGameofLifeClass, 1ProcStillLife)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string meshName1 = "1ProcPartStillLifeTest";
        std::string meshName2 = "1ProcFieldStillLifeTest";

        QuadGameofLifeMesh PartMesh(comm, 4, 4);
        QuadGameofLifeMesh FieldMesh(comm, 4, 4);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        stk::mesh::EntityIdVector elemIds = {6, 7, 10, 11};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        ASSERT_EQ(4u, PartGame.get_num_active_elements());
        ASSERT_EQ(4u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(10);
        FieldGame.run_game_of_life(10);

        ASSERT_EQ(4u, PartGame.get_num_active_elements());
        ASSERT_EQ(4u, FieldGame.get_num_active_elements());
    }
}
TEST(QuadGameofLifeClass, 4ProcStillLife)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string meshName1 = "4ProcPartStillLifeTest";
        std::string meshName2 = "4ProcFieldStillLifeTest";

        QuadGameofLifeMesh PartMesh(comm, 4, 4);
        QuadGameofLifeMesh FieldMesh(comm, 4, 4);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        stk::mesh::EntityIdVector elemIds = {6, 7, 10, 11};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        if (1 == procRank || 2 == procRank)
        {
            ASSERT_EQ(2u, PartGame.get_num_active_elements());
            ASSERT_EQ(2u, FieldGame.get_num_active_elements());
        }
        else
        {
            ASSERT_EQ(0u, PartGame.get_num_active_elements());
            ASSERT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(10);
        FieldGame.run_game_of_life(10);

        if (1 == procRank || 2 == procRank)
        {
            ASSERT_EQ(2u, PartGame.get_num_active_elements());
            ASSERT_EQ(2u, FieldGame.get_num_active_elements());
        }
        else
        {
            ASSERT_EQ(0u, PartGame.get_num_active_elements());
            ASSERT_EQ(0u, FieldGame.get_num_active_elements());
        }
    }
}
TEST(QuadGameofLifeClass, 1ProcOscillatorPeriod2)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string meshName1 = "1ProcPartOscillatorPeriod2";
        std::string meshName2 = "1ProcOscillatorFieldPeriod2";
        QuadGameofLifeMesh PartMesh(comm, 3, 3);
        QuadGameofLifeMesh FieldMesh(comm, 3, 3);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        stk::mesh::EntityIdVector elemIds = {4, 5, 6};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        EXPECT_EQ(3u, PartGame.get_num_active_elements());
        EXPECT_EQ(3u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector elemIds2 = {2, 5, 8};

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds2));
        EXPECT_EQ(3u, PartGame.get_num_active_elements());
        EXPECT_EQ(3u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector elemIds3 = {4, 5, 6};
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds3));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds3));
        EXPECT_EQ(3u, PartGame.get_num_active_elements());
        EXPECT_EQ(3u, FieldGame.get_num_active_elements());
    }

}
TEST(QuadGameofLifeClass, 4ProcOscillatorPeriod2)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);

    if (4 == numProcs)
    {
        std::string meshName1 = "4ProcPartOscillatorPeriod2";
        std::string meshName2 = "4ProcFieldOscillatorPeriod2";

        QuadGameofLifeMesh PartMesh(comm, 3, 3);
        QuadGameofLifeMesh FieldMesh(comm, 3, 3);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        stk::mesh::EntityIdVector elemIds = {4, 5, 6};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds));
        if (3 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector elemIds2 = {2, 5, 8};
        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds2));
        if (3 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (3 != procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds));
        if (3 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
    }
}
TEST(QuadGameofLifeClass, 1ProcOscillatorPeriod8)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string meshName1 = "1ProcPartOscillatorPeriod8";
        std::string meshName2 = "1ProcFieldOscillatorPeriod8";

        unsigned width = 10;
        unsigned rowsPerProc = 10;
        QuadGameofLifeMesh PartMesh(comm, width, rowsPerProc);
        QuadGameofLifeMesh FieldMesh(comm, width, rowsPerProc);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        stk::mesh::EntityIdVector elemIdsToActivate = {26, 27, 28, 36, 37, 38, 46, 47, 48,
                                                       53, 54, 55, 63, 64, 65, 73, 74, 75};
        PartGame.activate_these_ids(elemIdsToActivate);
        FieldGame.activate_these_ids(elemIdsToActivate);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIdsToActivate));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIdsToActivate));
        EXPECT_EQ(18u, PartGame.get_num_active_elements());
        EXPECT_EQ(18u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(4);
        FieldGame.run_game_of_life(4);

        stk::mesh::EntityIdVector exElemIds = {7, 17, 18, 26, 27, 29, 35, 38, 39, 40,
                                               44, 46, 48, 53, 55, 57, 61, 62, 63, 66,
                                               72, 74, 75, 83, 84, 94};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
        EXPECT_EQ(26u, PartGame.get_num_active_elements());
        EXPECT_EQ(26u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(4);
        FieldGame.run_game_of_life(4);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIdsToActivate));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIdsToActivate));
        EXPECT_EQ(18u, PartGame.get_num_active_elements());
        EXPECT_EQ(18u, FieldGame.get_num_active_elements());
    }
}
TEST(QuadGameofLifeClass, 4ProcOscillatorPeriod8)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procNum = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string meshName1 = "4ProcPartOscillatorPeriod8";
        std::string meshName2 = "4ProcFieldOscillatorPeriod8";
        unsigned width = 10;
        unsigned height = 10;

        QuadGameofLifeMesh PartMesh(comm, width, height);
        QuadGameofLifeMesh FieldMesh(comm, width, height);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        stk::mesh::EntityIdVector elemIdsToActivate = {26, 27, 28, 36, 37, 38, 46, 47, 48,
                                                       53, 54, 55, 63, 64, 65, 73, 74, 75};
        PartGame.activate_these_ids(elemIdsToActivate);
        FieldGame.activate_these_ids(elemIdsToActivate);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIdsToActivate));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIdsToActivate));
        if (0 != procNum)
        {
            EXPECT_EQ(6u, PartGame.get_num_active_elements());
            EXPECT_EQ(6u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(4);
        FieldGame.run_game_of_life(4);

        stk::mesh::EntityIdVector exElemIds = {7, 17, 18, 26, 27, 29, 35, 38, 39, 40,
                                               44, 46, 48, 53, 55, 57, 61, 62, 63, 66,
                                               72, 74, 75, 83, 84, 94};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
        if (0 == procNum)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (1 == procNum)
        {
            EXPECT_EQ(7u, PartGame.get_num_active_elements());
            EXPECT_EQ(7u, FieldGame.get_num_active_elements());
        }
        else if (2 == procNum)
        {
            EXPECT_EQ(6u, PartGame.get_num_active_elements());
            EXPECT_EQ(6u, FieldGame.get_num_active_elements());
        }
        else if (3 == procNum)
        {
            EXPECT_EQ(10u, PartGame.get_num_active_elements());
            EXPECT_EQ(10u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(4);
        FieldGame.run_game_of_life(4);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIdsToActivate));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIdsToActivate));
        if (0 != procNum)
        {
            EXPECT_EQ(6u, PartGame.get_num_active_elements());
            EXPECT_EQ(6u, FieldGame.get_num_active_elements());
        }
        else
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
    }
}
TEST(HexGameofLifeClass, 1ProcBasicTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string partMeshName = "1ProcPartBasic";
        std::string fieldMeshName = "1ProcFieldBasic";

        HexGameofLifeMesh PartMesh(comm, 3, 3, 4);
        HexGameofLifeMesh FieldMesh(comm, 3, 3, 4);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, fieldMeshName);

        EXPECT_EQ(36u, PartGame.get_num_elems_on_proc());
        EXPECT_EQ(36u, FieldGame.get_num_elems_on_proc());

        stk::mesh::EntityIdVector elemIdsToActivate = {1, 3, 7, 9, 14, 23, 28, 30, 34, 36};
        PartGame.activate_these_ids(elemIdsToActivate);
        FieldGame.activate_these_ids(elemIdsToActivate);
        EXPECT_EQ(10u, PartGame.get_num_active_elements());
        EXPECT_EQ(10u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds = {5, 14, 23, 32};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
        EXPECT_EQ(4u, PartGame.get_num_active_elements());
        EXPECT_EQ(4u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        EXPECT_EQ(0u, PartGame.get_num_active_elements());
        EXPECT_EQ(0u, FieldGame.get_num_active_elements());
    }
}
TEST(HexGameofLifeClass, 4ProcBasicTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string partMeshName = "4ProcPartBasic";
        std::string fieldMeshName = "4ProcFieldBasic";

        HexGameofLifeMesh PartMesh(comm, 3, 3, 4);
        HexGameofLifeMesh FieldMesh(comm, 3, 3, 4);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, fieldMeshName);

        EXPECT_EQ(9u, PartGame.get_num_elems_on_proc());
        EXPECT_EQ(9u, FieldGame.get_num_elems_on_proc());

        stk::mesh::EntityIdVector elemIdsToActivate = {1, 3, 7, 9, 14, 23, 28, 30, 34, 36};
        PartGame.activate_these_ids(elemIdsToActivate);
        FieldGame.activate_these_ids(elemIdsToActivate);
        if (0 == procRank || 3 == procRank)
        {
            EXPECT_EQ(4u, PartGame.get_num_active_elements());
            EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank || 2 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds = {5, 14, 23, 32};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
        EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        EXPECT_EQ(1u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        EXPECT_EQ(0u, PartGame.get_num_active_elements());
        EXPECT_EQ(0u, FieldGame.get_num_active_elements());
    }
}
TEST(HexGameofLifeClass, 1ProcOscillator)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string partMeshName = "1ProcPartOscillator";
        std::string fieldMeshName = "1ProcFieldOscillator";
        HexGameofLifeMesh PartMesh(comm, 3, 3, 4);
        HexGameofLifeMesh FieldMesh(comm, 3, 3, 4);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, fieldMeshName);


        EXPECT_EQ(36u, PartGame.get_num_elems_on_proc());
        EXPECT_EQ(36u, FieldGame.get_num_elems_on_proc());

        stk::mesh::EntityIdVector elemIds = {10, 12, 13, 15, 17, 19, 21, 22, 24, 26};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds));
        EXPECT_EQ(10u, PartGame.get_num_active_elements());
        EXPECT_EQ(10u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds1 = {5, 13, 15, 17, 22, 24, 26, 32};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds1));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds1));
        EXPECT_EQ(8u, PartGame.get_num_active_elements());
        EXPECT_EQ(8u, PartGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds2 = {11, 13, 15, 16, 18, 20, 22, 24, 25, 27};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds2));
        EXPECT_EQ(10u, PartGame.get_num_active_elements());
        EXPECT_EQ(10u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds3 = {5, 11, 13, 15, 20, 22, 24, 32};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds3));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds3));
        EXPECT_EQ(8u, PartGame.get_num_active_elements());
        EXPECT_EQ(8u, PartGame.get_num_active_elements());

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds));
        EXPECT_EQ(10u, PartGame.get_num_active_elements());
        EXPECT_EQ(10u, FieldGame.get_num_active_elements());
    }
}
TEST(HexGameofLifeClass, 4ProcOscillator)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string partMeshName = "4ProcPartOscillator";
        std::string fieldMeshName = "4ProcFieldOscillator";
        HexGameofLifeMesh PartMesh(comm, 3, 3, 4);
        HexGameofLifeMesh FieldMesh(comm, 3, 3, 4);

        PartGameofLife PartGame(PartMesh, partMeshName);
        FieldGameofLife FieldGame(FieldMesh, fieldMeshName);

        EXPECT_EQ(9u, PartGame.get_num_elems_on_proc());
        EXPECT_EQ(9u, FieldGame.get_num_elems_on_proc());

        stk::mesh::EntityIdVector elemIds = {10, 12, 13, 15, 17, 19, 21, 22, 24, 26};
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds));
        if (3 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else if (0 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds1 = {5, 13, 15, 17, 22, 24, 26, 32};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds1));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds1));
        if (3 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (0 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds2 = {11, 13, 15, 16, 18, 20, 22, 24, 25, 27};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds2));
        if (3 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else if (0 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        stk::mesh::EntityIdVector exElemIds3 = {5, 11, 13, 15, 20, 22, 24, 32};
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds3));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds3));
        if (3 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (0 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }

        PartGame.run_game_of_life(1);
        FieldGame.run_game_of_life(1);

        EXPECT_TRUE(PartGame.are_these_ids_active(elemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(elemIds));
        if (3 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else if (0 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
            EXPECT_EQ(0u, FieldGame.get_num_active_elements());
        }
    }
}
/*
 *
 *
 *
 */
TEST(PNG, 1ProcBordered)
{
    std::string fileName = "Boss.png";
    BorderedPNGProcessor png(fileName);
    png.commit_image_vector_to_pixel_vector();

    EXPECT_EQ(13u, png.get_image_width());
    EXPECT_EQ(16u, png.get_image_height());

}
TEST(PNG, 1ProcPaddingImages)
{
    std::string fileName = "Tiny.png";
    PNGProcessor png(fileName);
    png.commit_image_vector_to_pixel_vector();

    EXPECT_EQ(9u, png.get_image_width());
    EXPECT_EQ(10u, png.get_image_height());

    png.add_this_much_pixel_padding_to_bottom(2);
    EXPECT_EQ(12u, png.get_image_height());

    png.add_this_much_pixel_padding_to_top(2);
    EXPECT_EQ(14u, png.get_image_height());

    png.add_this_much_pixel_padding_to_left(2);
    EXPECT_EQ(11u, png.get_image_width());

    png.add_this_much_pixel_padding_to_right(2);
    EXPECT_EQ(13u, png.get_image_width());
}
TEST(PNGGameofLife, 1ProcCompressionTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string fileName = "Aircraftcarrier.png";
        std::string meshName = "1ProcCarrierCompression";

        BorderedPNGProcessor carrier(fileName);
        carrier.commit_image_vector_to_pixel_vector();

        stk::mesh::EntityIdVector elemIds;
        carrier.fill_id_vector_with_active_pixels(elemIds);

        unsigned width = carrier.get_image_width();
        unsigned height = carrier.get_image_height();

        QuadGameofLifeMesh Mesh(comm, width, height);

        FieldGameofLife FieldGame(Mesh, meshName);
        FieldGame.activate_these_ids(elemIds);

        stk::mesh::EntityIdVector exElemIds = {10, 11, 17, 14, 20,21};
        EXPECT_EQ(6u, FieldGame.get_num_active_elements());
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
    }
}
TEST(PNGGameofLife, 4ProcCompressionTest)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string fileName = "Aircraftcarrier.png";
        std::string meshName = "4ProcCarrierCompression";

        BorderedPNGProcessor carrier(fileName);
        carrier.commit_image_vector_to_pixel_vector();

        stk::mesh::EntityIdVector elemIds;
        carrier.fill_id_vector_with_active_pixels(elemIds);

        unsigned width = carrier.get_image_width();
        unsigned height = carrier.get_image_height();

        QuadGameofLifeMesh Mesh (comm, width, height);
        PartGameofLife PartGame(Mesh, meshName);
        PartGame.activate_these_ids(elemIds);

        stk::mesh::EntityIdVector exElemIds = {10, 11, 14, 17, 20, 21};

        if (0 == procRank)
        {
            EXPECT_EQ(0u, PartGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
        }
        else if (3 == procRank)
        {
            EXPECT_EQ(2u, PartGame.get_num_active_elements());
        }
    }
}
TEST(PNGGameofLife, 1ProcTiny)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string fileName = "Tiny.png";
        std::string meshName1 = "1ProcPartTiny";
        std::string meshName2 = "1ProcFieldTiny";

        PNGProcessor PNG(fileName);

        EXPECT_EQ(9u, PNG.get_image_width());
        EXPECT_EQ(10u, PNG.get_image_height());

        PNG.commit_image_vector_to_pixel_vector();

        EXPECT_EQ(9u, PNG.get_image_width());
        EXPECT_EQ(10u, PNG.get_image_height());

        unsigned width = PNG.get_image_width();
        unsigned height = PNG.get_image_height();

        QuadGameofLifeMesh PartMesh(comm, width, height);
        QuadGameofLifeMesh FieldMesh(comm, width, height);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName2);

        stk::mesh::EntityIdVector elemIds;
        PNG.fill_id_vector_with_active_pixels(elemIds);
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        stk::mesh::EntityIdVector exElemIds1 =
        {14, 22, 23, 24, 30, 31, 32, 33, 34, 38,
         39, 40, 41, 42, 43, 44, 48, 49, 50, 51,
         52, 58, 59, 60, 68};

        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds1));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds1));
        EXPECT_EQ(25u, PartGame.get_num_active_elements());
        EXPECT_EQ(25u, FieldGame.get_num_active_elements());

        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        stk::mesh::EntityIdVector exElemIds2 =
        {
         3, 7, 12, 16, 19, 20, 26, 27, 32,
         40, 42, 50, 55, 56, 62, 63, 66,
         70, 75, 79, 85, 86, 87
        };
        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds2));
        EXPECT_EQ(23u, PartGame.get_num_active_elements());
        EXPECT_EQ(23u, FieldGame.get_num_active_elements());
    }
}

enum PixelColor { RED = 2, GREEN, BLUE };

std::vector<Pixel> get_colored_pixels_by_color(SimpleColoredPng & image, enum PixelColor pixelColor)
{
    switch(pixelColor)
    {
        case RED:
            return image.get_red_color_coords();
        case GREEN:
            return image.get_green_color_coords();
        case BLUE:
            return image.get_blue_color_coords();
        default:
            break;
    }
    return {};
}


void create_nodeset_for_colored_pixels(stk::mesh::BulkData & bulk, SimpleColoredPng & image, enum PixelColor pixelColor)
{
    std::vector<Pixel> coloredPixels = get_colored_pixels_by_color(image, pixelColor);

    stk::mesh::EntityIdVector elementIds = image.get_elemIds_for_colored_pixels(coloredPixels);

    stk::mesh::EntityVector elementNodes;
    for(stk::mesh::EntityId elemId : elementIds)
    {
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
        elementNodes.insert(elementNodes.end(), bulk.begin_nodes(elem), bulk.end_nodes(elem));
    }

    bulk.modification_begin();
    std::string partName = "nodelist_" + std::to_string(pixelColor);
    stk::mesh::Part& nodesetPart = bulk.mesh_meta_data().declare_part(partName, stk::topology::NODE_RANK);
    for(stk::mesh::Entity node : elementNodes)
        bulk.change_entity_parts(node, stk::mesh::PartVector {&nodesetPart});
    bulk.modification_end();
}

TEST(TOSDTWD, quad_mesh_from_png)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string fileName = unitTestUtils::getOption("-i", "Tiny.png");
        SimpleColoredPng image(fileName);
        unsigned width = image.get_image_width();
        unsigned height = image.get_image_height();

        QuadGameofLifeMesh FieldMesh(comm, width, height);
        FieldGameofLife FieldGame(FieldMesh, "junk");

        stk::mesh::EntityIdVector elemIds;
        image.fill_id_vector_with_active_pixels(elemIds);
        FieldGame.activate_these_ids(elemIds);

        stk::mesh::BulkData &bulk = FieldMesh.bulk_data();
        create_nodeset_for_colored_pixels(bulk, image, RED);
        create_nodeset_for_colored_pixels(bulk, image, GREEN);
        create_nodeset_for_colored_pixels(bulk, image, BLUE);

        FieldGame.write_mesh();
    }
}

TEST(TOSDTWD, hex_mesh_from_png)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        std::string fileName = unitTestUtils::getOption("-i", "Tiny.png");
        SimpleColoredPng image(fileName);
        unsigned width = image.get_image_width();
        unsigned height = image.get_image_height();

        HexGameofLifeMesh FieldMesh(comm, width, height, 1);
        FieldGameofLife FieldGame(FieldMesh, "junk");

        stk::mesh::EntityIdVector elemIds;
        image.fill_id_vector_with_active_pixels(elemIds);
        FieldGame.activate_these_ids(elemIds);

        stk::mesh::BulkData &bulk = FieldMesh.bulk_data();
        create_nodeset_for_colored_pixels(bulk, image, RED);
        create_nodeset_for_colored_pixels(bulk, image, GREEN);
        create_nodeset_for_colored_pixels(bulk, image, BLUE);

        FieldGame.write_mesh();
    }
}


TEST(PNGGameofLife, 4ProcTiny)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        std::string fileName = "Tiny.png";
        std::string meshName1 = "4ProcPartTiny";
        std::string meshName2 = "4ProcFieldTiny";

        PNGProcessor PNG(fileName);

        EXPECT_EQ(9u, PNG.get_image_width());
        EXPECT_EQ(10u, PNG.get_image_height());

        PNG.commit_image_vector_to_pixel_vector();

        EXPECT_EQ(9u, PNG.get_image_width());
        EXPECT_EQ(10u, PNG.get_image_height());

        unsigned width = PNG.get_image_width();
        unsigned height = PNG.get_image_height();

        QuadGameofLifeMesh PartMesh(comm, width, height);
        QuadGameofLifeMesh FieldMesh(comm, width, height);

        PartGameofLife PartGame(PartMesh, meshName1);
        FieldGameofLife FieldGame(FieldMesh, meshName1);

        stk::mesh::EntityIdVector elemIds;
        PNG.fill_id_vector_with_active_pixels(elemIds);
        PartGame.activate_these_ids(elemIds);
        FieldGame.activate_these_ids(elemIds);

        stk::mesh::EntityIdVector exElemIds = {14, 22, 23, 24, 30, 31, 32, 33, 34, 38,
                                               39, 40, 41, 42, 43, 44, 48, 49, 50, 51,
                                               52, 58, 59, 50, 58};

        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds));
        if (0 == procRank)
        {
            EXPECT_EQ(1u, PartGame.get_num_active_elements());
            EXPECT_EQ(1u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(8u, PartGame.get_num_active_elements());
            EXPECT_EQ(8u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(12u, PartGame.get_num_active_elements());
            EXPECT_EQ(12u, FieldGame.get_num_active_elements());
        }
        else if (3 == procRank)
        {
            EXPECT_EQ(4u, PartGame.get_num_active_elements());
            EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        }


        PartGame.run_game_of_life(5);
        FieldGame.run_game_of_life(5);

        stk::mesh::EntityIdVector exElemIds2 = {3, 7, 12, 16, 19, 20, 26, 27, 32, 40, 42, 50,
                                                55, 56, 62, 63, 66, 70, 75, 79};

        EXPECT_TRUE(PartGame.are_these_ids_active(exElemIds2));
        EXPECT_TRUE(FieldGame.are_these_ids_active(exElemIds2));
        if (0 == procRank)
        {
            EXPECT_EQ(4u, PartGame.get_num_active_elements());
            EXPECT_EQ(4u, FieldGame.get_num_active_elements());
        }
        else if (1 == procRank)
        {
            EXPECT_EQ(5u, PartGame.get_num_active_elements());
            EXPECT_EQ(5u, FieldGame.get_num_active_elements());
        }
        else if (2 == procRank)
        {
            EXPECT_EQ(3u, PartGame.get_num_active_elements());
            EXPECT_EQ(3u, FieldGame.get_num_active_elements());
        }
        else if (3 == procRank)
        {
            EXPECT_EQ(11u, PartGame.get_num_active_elements());
            EXPECT_EQ(11u, FieldGame.get_num_active_elements());
        }
    }
}
TEST(NewNoGhostGame, 1ProcGeneralStuff)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexGameofLifeMesh Mesh(comm, 1, 1, 1, stk::mesh::BulkData::NO_AUTO_AURA);
        NoGhostGameofLife PacMan(Mesh, "1ProcGeneralStuff");

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(PacMan.bulk_data(), stk::topology::ELEM_RANK, elements);
        stk::mesh::Entity elem1 = PacMan.element_with_id(1);

        EXPECT_EQ(1u, elements.size());
        EXPECT_TRUE(PacMan.is_valid_entity(elem1));
        EXPECT_EQ(1u, PacMan.num_procs());
        EXPECT_EQ(1u, PacMan.num_elems_on_proc());
    }
}
TEST(NewNoGhostGame, 4ProcGeneralStuff)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexGameofLifeMesh Mesh(comm, 1, 1, 4, stk::mesh::BulkData::NO_AUTO_AURA);
        NoGhostGameofLife PacMan(Mesh, "4ProcGeneralStuff");

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(PacMan.bulk_data(), stk::topology::ELEM_RANK, elements);

        if (0 == procRank)
        {
            stk::mesh::Entity elem = PacMan.element_with_id(1);
            EXPECT_TRUE(PacMan.is_valid_entity(elem));
        }
        if (1 == procRank)
        {
            stk::mesh::Entity elem = PacMan.element_with_id(2);
            EXPECT_TRUE(PacMan.is_valid_entity(elem));
        }
        if (2 == procRank)
        {
            stk::mesh::Entity elem = PacMan.element_with_id(3);
            EXPECT_TRUE(PacMan.is_valid_entity(elem));
        }
        if (3 == procRank)
        {
            stk::mesh::Entity elem = PacMan.element_with_id(4);
            EXPECT_TRUE(PacMan.is_valid_entity(elem));
        }

        EXPECT_EQ(1u, elements.size());
        EXPECT_EQ(4u, PacMan.num_procs());
        EXPECT_EQ(1u, PacMan.num_elems_on_proc());
    }
}
TEST(NewNoGhostGame, 1ProcLessGeneralStuff)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexGameofLifeMesh Mesh(comm, 2, 2, 2, stk::mesh::BulkData::NO_AUTO_AURA);
        NoGhostGameofLife PacMan(Mesh, "1ProcLessGeneralStuff");

        for (unsigned id = 1; id <= 8; id++)
        {
            stk::mesh::Entity elem = PacMan.element_with_id(id);
            EXPECT_EQ(7u, PacMan.num_neighbors(elem));
        }

        EXPECT_EQ(0u, PacMan.num_active_elems());

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4, 5};
        PacMan.activate_these_ids(elemIds);

        EXPECT_EQ(5u, PacMan.num_active_elems());
    }
}
TEST(NewNoGhostGame, 4ProcLessGeneralStuff)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexGameofLifeMesh Mesh(comm, 1, 1, 8, stk::mesh::BulkData::NO_AUTO_AURA);
        NoGhostGameofLife PacMan(Mesh, "4ProcLessGeneralStuff4");

        if (0 == procRank)
        {
            stk::mesh::Entity elem1 = PacMan.element_with_id(1);
            stk::mesh::Entity elem2 = PacMan.element_with_id(2);
            EXPECT_EQ(1u, PacMan.num_neighbors(elem1));
            EXPECT_EQ(2u, PacMan.num_neighbors(elem2));
        }
        if (1 == procRank)
        {
            stk::mesh::Entity elem1 = PacMan.element_with_id(3);
            stk::mesh::Entity elem2 = PacMan.element_with_id(4);
            EXPECT_EQ(2u, PacMan.num_neighbors(elem1));
            EXPECT_EQ(2u, PacMan.num_neighbors(elem2));
        }
        if (2 == procRank)
        {
            stk::mesh::Entity elem1 = PacMan.element_with_id(5);
            stk::mesh::Entity elem2 = PacMan.element_with_id(6);
            EXPECT_EQ(2u, PacMan.num_neighbors(elem1));
            EXPECT_EQ(2u, PacMan.num_neighbors(elem2));
        }
        if (3 == procRank)
        {
            stk::mesh::Entity elem1 = PacMan.element_with_id(7);
            stk::mesh::Entity elem2 = PacMan.element_with_id(8);
            EXPECT_EQ(2u, PacMan.num_neighbors(elem1));
            EXPECT_EQ(1u, PacMan.num_neighbors(elem2));
        }
        EXPECT_EQ(0u, PacMan.num_active_elems());

        stk::mesh::EntityIdVector elemIds = {1, 3, 5, 7};
        PacMan.activate_these_ids(elemIds);

        EXPECT_EQ(1u, PacMan.num_active_elems());
    }
}
TEST(NewNoGhostGame, 1ProcNeighbors)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);

    if (1 == numProcs)
    {
        HexGameofLifeMesh Mesh(MPI_COMM_WORLD, 2, 2, 2);
        NoGhostGameofLife Game(Mesh, "1ProcNeighbors");

        stk::mesh::Entity elem1 = Game.element_with_id(1);
        stk::mesh::Entity elem8 = Game.element_with_id(8);

        EXPECT_EQ(7u, Game.num_neighbors(elem1));
        EXPECT_EQ(7u, Game.num_neighbors(elem8));

        EXPECT_EQ(0u, Game.num_active_neighbors(elem1));
        EXPECT_EQ(0u, Game.num_active_neighbors(elem8));

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4};
        Game.activate_these_ids(elemIds);

        EXPECT_EQ(3u, Game.num_active_neighbors(elem1));
        EXPECT_EQ(4u, Game.num_active_neighbors(elem8));
    }
}
TEST(NewNoGhostGame, 4ProcNeighbors)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procNum = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexGameofLifeMesh Mesh(MPI_COMM_WORLD, 2, 2, 4, stk::mesh::BulkData::NO_AUTO_AURA);
        NoGhostGameofLife Game(Mesh, "4ProcNeighbors");

        stk::mesh::Entity elem;
        unsigned expectedNeighbors;
        if (0 == procNum)
        {
            elem = Game.element_with_id(1);
            expectedNeighbors = 3u;
        }
        else if (1 == procNum)
        {
            elem = Game.element_with_id(5);
            expectedNeighbors = 8u;
        }
        else if (2 == procNum)
        {
            elem = Game.element_with_id(9);
            expectedNeighbors = 4u;
        }
        else if (3 == procNum)
        {
            elem = Game.element_with_id(13);
            expectedNeighbors = 5u;
        }

        EXPECT_EQ(0u, Game.num_active_neighbors(elem));

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4, 9, 10, 11, 12, 14};
        Game.activate_these_ids(elemIds);

        EXPECT_EQ(expectedNeighbors, Game.num_active_neighbors(elem));
    }
}
TEST(NewNoGhostGame, 1ProcRunGame)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        HexGameofLifeMesh Mesh(comm, 2, 2, 2);

        NoGhostGameofLife Game(Mesh, "1ProcRunGame");

        stk::mesh::EntityIdVector elemIds = {1, 2, 3, 4, 5};
        Game.activate_these_ids(elemIds);
        EXPECT_EQ(5u, Game.num_active_elems());

        Game.run_game_of_life(1);
        EXPECT_EQ(8u, Game.num_active_elems());

        Game.run_game_of_life(1);
        EXPECT_EQ(0u, Game.num_active_elems());

        Game.run_game_of_life(1);
        EXPECT_EQ(0u, Game.num_active_elems());
    }
}
TEST(NewNoGhostGame, 4ProcRunGame)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procNum = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        HexGameofLifeMesh Mesh(comm, 3, 3, 4, stk::mesh::BulkData::NO_AUTO_AURA);

        NoGhostGameofLife Game(Mesh, "4ProcRunGame");

        stk::mesh::EntityIdVector elemIds = {10, 12, 13, 15, 17, 19, 21, 22, 24, 26};
        Game.activate_these_ids(elemIds);

        unsigned expectedActive;
        if (0 == procNum)
            expectedActive = 0;
        else if (1 == procNum)
            expectedActive = 5;
        else if (2 == procNum)
            expectedActive = 5;
        else if (3 == procNum)
            expectedActive = 0;

        EXPECT_EQ(expectedActive, Game.num_active_elems());

        Game.run_game_of_life(1);

        if (0 == procNum)
            expectedActive = 1;
        else if (1 == procNum)
            expectedActive = 3;
        else if (2 == procNum)
            expectedActive = 3;
        else if (3 == procNum)
            expectedActive = 1;

        EXPECT_EQ(expectedActive, Game.num_active_elems());

        Game.run_game_of_life(1);

        if (0 == procNum)
            expectedActive = 0;
        else if (1 == procNum)
            expectedActive = 5;
        else if (2 == procNum)
            expectedActive = 5;
        else if (3 == procNum)
            expectedActive = 0;

        EXPECT_EQ(expectedActive, Game.num_active_elems());

        Game.run_game_of_life(1);

        if (0 == procNum)
            expectedActive = 1;
        else if (1 == procNum)
            expectedActive = 3;
        else if (2 == procNum)
            expectedActive = 3;
        else if (3 == procNum)
            expectedActive = 1;

        EXPECT_EQ(expectedActive, Game.num_active_elems());
    }
}
TEST(NewNoGhostGame, 1ProcQuad)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
        QuadGameofLifeMesh Mesh(comm, 8, 8);

        NoGhostGameofLife Game(Mesh, "1ProcQuad");
        stk::mesh::EntityIdVector elemIds = {41, 42, 43, 51, 58};
        Game.activate_these_ids(elemIds);

        EXPECT_EQ(5u, Game.num_active_elems());

        stk::mesh::Entity elem1 = Game.element_with_id(50);
        stk::mesh::Entity elem2 = Game.element_with_id(41);

        EXPECT_EQ(8u, Game.num_neighbors(elem1));
        EXPECT_EQ(5u, Game.num_active_neighbors(elem1));
        EXPECT_EQ(5u, Game.num_neighbors(elem2));
        EXPECT_EQ(1u, Game.num_active_neighbors(elem2));

        Game.run_game_of_life(2);

        EXPECT_EQ(5u, Game.num_active_elems());

        EXPECT_EQ(8u, Game.num_neighbors(elem1));
        EXPECT_EQ(3u, Game.num_active_neighbors(elem1));
        EXPECT_EQ(5u, Game.num_neighbors(elem2));
        EXPECT_EQ(1u, Game.num_active_neighbors(elem2));

        Game.run_game_of_life(2);

        stk::mesh::Entity elem3 = Game.element_with_id(36);
        stk::mesh::Entity elem4 = Game.element_with_id(28);

        EXPECT_EQ(8u, Game.num_neighbors(elem3));
        EXPECT_EQ(2u, Game.num_active_neighbors(elem3));
        EXPECT_EQ(8u, Game.num_neighbors(elem4));
        EXPECT_EQ(2u, Game.num_active_neighbors(elem4));

        Game.run_game_of_life(2);


        EXPECT_EQ(8u, Game.num_neighbors(elem3));
        EXPECT_EQ(3u, Game.num_active_neighbors(elem3));
        EXPECT_EQ(8u, Game.num_neighbors(elem4));
        EXPECT_EQ(2u, Game.num_active_neighbors(elem4));
    }
}
TEST(NewNoGhostGame, 4ProcQuad)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
        QuadGameofLifeMesh Mesh(comm, 8, 8, stk::mesh::BulkData::NO_AUTO_AURA);

        NoGhostGameofLife Game(Mesh, "4ProcQuad");
        stk::mesh::EntityIdVector elemIds = {41, 42, 43, 51, 58};
        Game.activate_these_ids(elemIds);

        stk::mesh::Entity elem;
        unsigned expectedActive;
        unsigned expectedActiveNeighbors;
        if (0 == procRank)
        {
            expectedActive = 0;
            elem = Game.element_with_id(9);
            expectedActiveNeighbors = 0;
        }
        else if (1 == procRank)
        {
            expectedActive = 0;
            elem = Game.element_with_id(25);
            expectedActiveNeighbors = 0;
        }
        else if (2 == procRank)
        {
            expectedActive = 3;
            elem = Game.element_with_id(41);
            expectedActiveNeighbors = 1;
        }
        else if (3 == procRank)
        {
            expectedActive = 2;
            elem = Game.element_with_id(57);
            expectedActiveNeighbors = 1;
        }

        EXPECT_EQ(expectedActive, Game.num_active_elems());
        EXPECT_EQ(expectedActiveNeighbors, Game.num_active_neighbors(elem));

        Game.run_game_of_life(2);

        if (0 == procRank)
        {
            expectedActive = 0;
            expectedActiveNeighbors = 0;
        }
        else if (1 == procRank)
        {
            expectedActive = 0;
            expectedActiveNeighbors = 1;
        }
        else if (2 == procRank)
        {
            expectedActive = 4;
            expectedActiveNeighbors = 1;
        }
        else if (3 == procRank)
        {
            expectedActive =1;
            expectedActiveNeighbors = 0;
        }

        EXPECT_EQ(expectedActive, Game.num_active_elems());
        EXPECT_EQ(expectedActiveNeighbors, Game.num_active_neighbors(elem));

        Game.run_game_of_life(2);

        if (0 == procRank)
        {
            expectedActive = 0;
            expectedActiveNeighbors = 0;
            elem = Game.element_with_id(11);
        }
        else if (1 == procRank)
        {
            expectedActive = 0;
            expectedActiveNeighbors = 3;
            elem = Game.element_with_id(27);
        }
        else if (2 == procRank)
        {
            expectedActive = 4;
            expectedActiveNeighbors = 5;
            elem = Game.element_with_id(43);
        }
        else if (3 == procRank)
        {
            expectedActive = 1;
            expectedActiveNeighbors = 1;
            elem = Game.element_with_id(59);
        }

        EXPECT_EQ(expectedActive, Game.num_active_elems());
        EXPECT_EQ(expectedActiveNeighbors, Game.num_active_neighbors(elem));

        Game.run_game_of_life(2);

        if (0 == procRank)
        {
            expectedActive = 0;
            expectedActiveNeighbors = 0;
        }
        else if (1 == procRank)
        {
            expectedActive = 2;
            expectedActiveNeighbors = 3;
        }
        else if (2 == procRank)
        {
            expectedActive = 3;
            expectedActiveNeighbors = 3;
        }
        else if (3 == procRank)
        {
            expectedActive = 0;
            expectedActiveNeighbors = 0;
        }

        EXPECT_EQ(expectedActive, Game.num_active_elems());
        EXPECT_EQ(expectedActiveNeighbors, Game.num_active_neighbors(elem));
    }
}
TEST(NewNoGhostGame, 1ProcTri)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (1 == numProcs)
    {
       TriGameofLifeMesh Mesh(comm, 4, 4, stk::mesh::BulkData::NO_AUTO_AURA);

       NoGhostGameofLife Game(Mesh, "1ProcTri");
       stk::mesh::EntityIdVector elemIds  = {1, 2, 3, 4, 5, 6, 7,
                                             8, 25, 26, 27, 28, 29, 30, 31, 32};
       Game.activate_these_ids(elemIds);

       EXPECT_EQ(16u, Game.num_active_elems());

       stk::mesh::Entity elem1 = Game.element_with_id(1);
       stk::mesh::Entity elem2 = Game.element_with_id(8);
       stk::mesh::Entity elem3 = Game.element_with_id(25);
       stk::mesh::Entity elem4 = Game.element_with_id(32);

       EXPECT_EQ(3u, Game.num_neighbors(elem1));
       EXPECT_EQ(6u, Game.num_neighbors(elem2));
       EXPECT_EQ(6u, Game.num_neighbors(elem3));
       EXPECT_EQ(3u, Game.num_neighbors(elem4));

       EXPECT_EQ(2u, Game.num_active_neighbors(elem1));
       EXPECT_EQ(2u, Game.num_active_neighbors(elem2));
       EXPECT_EQ(2u, Game.num_active_neighbors(elem3));
       EXPECT_EQ(2u, Game.num_active_neighbors(elem4));

       Game.run_game_of_life(2);

       EXPECT_EQ(14u, Game.num_active_elems());

       EXPECT_EQ(3u, Game.num_neighbors(elem1));
       EXPECT_EQ(6u, Game.num_neighbors(elem2));
       EXPECT_EQ(6u, Game.num_neighbors(elem3));
       EXPECT_EQ(3u, Game.num_neighbors(elem4));

       EXPECT_EQ(3u, Game.num_active_neighbors(elem1));
       EXPECT_EQ(2u, Game.num_active_neighbors(elem2));
       EXPECT_EQ(2u, Game.num_active_neighbors(elem3));
       EXPECT_EQ(3u, Game.num_active_neighbors(elem4));

       Game.run_game_of_life(2);

       EXPECT_EQ(6u, Game.num_active_elems());

       EXPECT_EQ(0u, Game.num_active_neighbors(elem1));
       EXPECT_EQ(1u, Game.num_active_neighbors(elem2));
       EXPECT_EQ(1u, Game.num_active_neighbors(elem3));
       EXPECT_EQ(0u, Game.num_active_neighbors(elem4));
    }
}
TEST(NewNoGhostGame, 4ProcTri)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (4 == numProcs)
    {
       TriGameofLifeMesh Mesh(comm, 4, 4, stk::mesh::BulkData::NO_AUTO_AURA);

       NoGhostGameofLife Game(Mesh, "4ProcTri");
       stk::mesh::EntityIdVector elemIds  = {1, 2, 3, 4, 5, 6, 7, 8, 25, 26,
                                             27, 28, 29, 30, 31, 32};
       Game.activate_these_ids(elemIds);

       unsigned expectedActive;
       stk::mesh::Entity elem;
       unsigned expectedNeighbors;
       unsigned expectedActiveNeighbors;

       if (0 == procRank)
       {
           expectedActive = 8;
           elem = Game.element_with_id(4);
           expectedNeighbors = 9;
           expectedActiveNeighbors = 4;
       }
       else if (1 == procRank)
       {
           expectedActive = 0;
           elem = Game.element_with_id(13);
           expectedNeighbors = 12;
           expectedActiveNeighbors = 5;
       }
       else if (2 == procRank)
       {
           expectedActive = 0;
           elem = Game.element_with_id(20);
           expectedNeighbors = 12;
           expectedActiveNeighbors = 5;

       }
       else if (3 == procRank)
       {
           expectedActive = 8;
           elem = Game.element_with_id(29);
           expectedNeighbors = 9;
           expectedActiveNeighbors = 4;
       }

       EXPECT_EQ(expectedActive, Game.num_active_elems());
       EXPECT_EQ(expectedNeighbors, Game.num_neighbors(elem));
       EXPECT_EQ(expectedActiveNeighbors, Game.num_active_neighbors(elem));

       Game.run_game_of_life(2);

       if (0 == procRank)
       {
           expectedActive = 5;
           expectedActiveNeighbors = 3;
       }
       else if (1 == procRank)
       {
           expectedActive = 2;
           expectedActiveNeighbors = 5;
       }
       else if (2 == procRank)
       {
           expectedActive = 2;
           expectedActiveNeighbors = 5;

       }
       else if (3 == procRank)
       {
           expectedActive = 5;
           expectedActiveNeighbors = 3;
       }

       EXPECT_EQ(expectedActive, Game.num_active_elems());
       EXPECT_EQ(expectedNeighbors, Game.num_neighbors(elem));
       EXPECT_EQ(expectedActiveNeighbors, Game.num_active_neighbors(elem));

       Game.run_game_of_life(2);

       if (0 == procRank)
       {
           expectedActive = 1;
           expectedActiveNeighbors = 1;
       }
       else if (1 == procRank)
       {
           expectedActive = 2;
           expectedActiveNeighbors = 4;
       }
       else if (2 == procRank)
       {
           expectedActive = 2;
           expectedActiveNeighbors = 4;

       }
       else if (3 == procRank)
       {
           expectedActive = 1;
           expectedActiveNeighbors = 1;
       }
    }
}
}




