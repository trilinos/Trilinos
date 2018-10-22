#include "gtest/gtest.h"
#include "mpi.h"
#include <string>
#include <sstream>
#include <vector>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_unit_test_utils/MeshFixture.hpp>
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

class TestTextMesh : public stk::unit_test_util::MeshFixture
{
protected:
    TestTextMesh()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    }
};

TEST_F(TestTextMesh, singleHex)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    }
}

void verify_coordinates(const std::vector<double> & coordinates, stk::mesh::BulkData & bulk, unsigned spatialDim)
{
    stk::mesh::EntityVector nodes;
    bulk.get_entities(stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part(), nodes);
    stk::unit_test_util::CoordinatesField & coordsField = static_cast<stk::unit_test_util::CoordinatesField&>(*bulk.mesh_meta_data().get_field(stk::topology::NODE_RANK, "coordinates"));
    for(size_t nodeIndex=0; nodeIndex < nodes.size(); nodeIndex++)
    {
       double * nodalCoords = stk::mesh::field_data(coordsField, nodes[nodeIndex]);
       for(unsigned coordIndex=0; coordIndex < spatialDim; coordIndex++)
           EXPECT_NEAR(coordinates[nodeIndex*spatialDim+coordIndex], nodalCoords[coordIndex], 1.0e-9);
    }
}

TEST_F(TestTextMesh, singleHexWithCoordinates)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    if (get_bulk().parallel_size() == 1)
    {
       std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };
       stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());
       verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
       verify_coordinates(coordinates, get_bulk(), 3);

    }
}

struct ElementInfo
{
    std::string blockName;
    stk::mesh::EntityId id;
};

void verify_part_membership(std::vector<ElementInfo> golds, const stk::mesh::MetaData &meta, const stk::mesh::BulkData &bulk)
{
    for(const ElementInfo gold : golds)
    {
        stk::mesh::Part *blockPart = meta.get_part(gold.blockName);
        ASSERT_TRUE(blockPart != nullptr);
        stk::mesh::EntityVector elems;
        stk::mesh::get_selected_entities(*blockPart, bulk.buckets(stk::topology::ELEM_RANK), elems);
        ASSERT_EQ(1u, elems.size());
        EXPECT_EQ(gold.id, bulk.identifier(elems[0]));
    }
}

TEST_F(TestTextMesh, twoHexDisconnectedWithCoordinatesAndParts)
{
    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n\
                            0,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
    if (get_bulk().parallel_size() == 1)
    {
       std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                           0,0,1, 1,0,1, 1,1,1, 0,1,1, 0,0,2, 1,0,2, 1,1,2, 0,1,2 };
       stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());
       verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
       verify_coordinates(coordinates, get_bulk(), 3);

       verify_part_membership({{"block_1", 1u}, {"block_2", 2u}}, get_meta(), get_bulk());
    }
}

TEST_F(TestTextMesh, mixedSpatialDim)
{
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         0,2,QUAD_4_2D,5,6,7,8";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TestTextMesh, singlHexWithSpaces)
{
    std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    }
}

TEST_F(TestTextMesh, singlHexWithLowerCase)
{
    std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    }
}

TEST_F(TestTextMesh, tooFewNodes)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TestTextMesh, tooManyNodes)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData)
{
    std::string meshDesc = "0,1,";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}

TEST_F(TestTextMesh, invalidTopology)
{
    std::string meshDesc = "0,1,invalid,1";
    EXPECT_THROW(stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk()), std::logic_error);
}


TEST_F(TestTextMesh, twoHexesSerial)
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

TEST_F(TestTextMesh, twoTet10Serial)
{
    std::string meshDesc =
        "0,1,TET_10,1,2,3,4,5,6,7,8,9,10\n\
         0,2,TET_10,2,11,3,4,12,13,6,9,14,10";
    if (get_bulk().parallel_size() == 1)
    {
  //                                         1       2      3        4          5          6         7           8 
        std::vector<double> coordinates = { 0,0,0, 1,0,0, 0.5,1,0, 0.5,0.5,1, 0.5,0,0, 0.75,0.5,0, 0.25,0.5,0, 0.25,0.25,0.5,
//                                           9              10            11         12           13        14
                                           0.75,0.25,0.5, 0.5,0.75,0.5, 1.5,0.5,0, 1.25,0.25,0, 1,0.75,0, 1,0.5,0.5 };
        stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());
        EXPECT_EQ(1,get_bulk().parallel_size());
        verify_single_element(get_bulk(), 1u, stk::topology::TET_10, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8,9,10});
        verify_single_element(get_bulk(), 2u, stk::topology::TET_10, stk::mesh::EntityIdVector{2,11,3,4,12,13,6,9,14,10});
        stk::io::write_mesh("twoTet10s.g", get_bulk());
    }
}

TEST_F(TestTextMesh, twoHexesParallel)
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

class TestTextMesh2d : public stk::unit_test_util::MeshFixture
{
protected:
    TestTextMesh2d() : stk::unit_test_util::MeshFixture(2)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    }
};

TEST_F(TestTextMesh2d, singleQuad)
{
    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
    if (get_bulk().parallel_size() == 1)
    {
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
        verify_single_element(get_bulk(), 1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
    }
}

TEST_F(TestTextMesh2d, threeQuadsWithCoordinates)
{
    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n\
                            0,2,QUAD_4_2D,2,3,5,6\n\
                            0,3,QUAD_4_2D,5,7,8,6";
    if (get_bulk().parallel_size() == 1)
    {
        std::vector<double> coordinates = { 0,0, 1,0, 1,1, 0,1, 2,0, 2,1, 3,0, 3,1, };
        stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

        std::vector<size_t> counts;
        stk::mesh::count_entities(get_meta().universal_part(), get_bulk(), counts);
        EXPECT_EQ(3u, counts[stk::topology::ELEM_RANK]);

        verify_coordinates(coordinates, get_bulk(), 2);
    }
}

TEST_F(TestTextMesh2d, DISABLED_twoQuadOneShellParallel)
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
