#include "gtest/gtest.h"
#include <mpi.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <string>
#include <sstream>
#include <vector>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>

namespace {

void verify_nodes(const stk::mesh::BulkData& bulk, stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds)
{
    stk::mesh::EntityVector nodes(bulk.begin_nodes(element), bulk.end_nodes(element));
    stk::mesh::EntityIdVector nodeIds(nodes.size());
    for(size_t i = 0; i < nodes.size(); ++i)
    {
        nodeIds[i] = bulk.identifier(nodes[i]);
    }
    EXPECT_EQ(goldNodeIds, nodeIds);
}

void verify_shared_nodes(const stk::mesh::BulkData& bulk, const stk::mesh::EntityIdVector& nodeIds , int sharingProc)
{
    for (stk::mesh::EntityId nodeId : nodeIds)
        EXPECT_TRUE(bulk.in_shared(stk::mesh::EntityKey(stk::topology::NODE_RANK,nodeId),sharingProc));
}

void verify_single_element(const stk::mesh::BulkData& bulk, stk::mesh::EntityId elemId, stk::topology topology, const stk::mesh::EntityIdVector& nodeIds)
{
    stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
    EXPECT_TRUE(bulk.is_valid(element));
    EXPECT_EQ(topology, bulk.bucket(element).topology());
    verify_nodes(bulk, element, nodeIds);
}

class TextMeshFixture : public stk::unit_test_util::MeshFixture
{
protected:
    TextMeshFixture()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    }
};

TEST_F(TextMeshFixture, singlHex)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    }
}

TEST_F(TextMeshFixture, singleQuad)
{
    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
    }
}

TEST_F(TextMeshFixture, mixedSpatialDim)
{
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         0,2,QUAD_4_2D,5,6,7,8";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TextMeshFixture, singlHexWithSpaces)
{
    std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    }
}

TEST_F(TextMeshFixture, singlHexWithLowerCase)
{
    std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    }
}

TEST_F(TextMeshFixture, tooFewNodes)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TextMeshFixture, tooManyNodes)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TextMeshFixture, tooLittleData)
{
    std::string meshDesc = "0,1,";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TextMeshFixture, invalidTopology)
{
    std::string meshDesc = "0,1,invalid,1";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}


TEST_F(TextMeshFixture, twoHexesSerial)
{
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         0,2,HEX_8,5,6,7,8,9,10,11,12";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        EXPECT_EQ(1,get_bulk().parallel_size());
        verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
        verify_single_element(get_bulk(), 2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
    }
}

TEST_F(TextMeshFixture, twoHexesParallel)
{
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         1,2,HEX_8,5,6,7,8,9,10,11,12";
    if (get_bulk().parallel_size() == 2)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        EXPECT_EQ(2,get_bulk().parallel_size());
        if (get_bulk().parallel_rank() == 0)
        {
            verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
            verify_shared_nodes(get_bulk(), {5,6,7,8}, 1);
        }
        else
        {
            verify_single_element(get_bulk(), 2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
            verify_shared_nodes(get_bulk(), {5,6,7,8}, 0);
        }
    }
}

TEST_F(TextMeshFixture, DISABLED_twoQuadOneShellParallel)
{
    std::string meshDesc =
        "0,1,QUAD_4_2D,1,2,3,4\n\
         1,2,QUAD_4_2D,3,4,5,6\n\
         2,3,SHELL_LINE_2,3,4";
    if (get_bulk().parallel_size() == 3)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        EXPECT_EQ(3,get_bulk().parallel_size());
        if (get_bulk().parallel_rank() == 0)
        {
            verify_single_element(get_bulk(), 1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
            verify_shared_nodes(get_bulk(), {3,4}, 1);
            verify_shared_nodes(get_bulk(), {3,4}, 2);
        }
        else if (get_bulk().parallel_rank() == 1)
        {
            verify_single_element(get_bulk(), 2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{3,4,5,6});
            verify_shared_nodes(get_bulk(), {3,4}, 0);
            verify_shared_nodes(get_bulk(), {3,4}, 2);
        }
        else
        {
            verify_single_element(get_bulk(), 3u, stk::topology::SHELL_LINE_2, stk::mesh::EntityIdVector{3,4});
            verify_shared_nodes(get_bulk(), {3,4}, 0);
            verify_shared_nodes(get_bulk(), {3,4}, 1);
        }
    }
}

} // namespace
