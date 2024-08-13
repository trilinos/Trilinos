#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc

#include <memory>
#include <vector>
#include <string>
#include <algorithm>

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_unit_test_utils/BuildMesh.hpp"
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include "stk_mesh/base/Comm.hpp"
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SidesetUpdater.hpp>
#include <stk_mesh/base/PolarityUtil.hpp>

#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp"

class TestTextMesh : public stk::unit_test_util::MeshFixture
{
protected:
  TestTextMesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void inititalize_2D_mesh()
  {
    reset_mesh();
    unsigned spatialDim = 2u;
    set_spatial_dimension(spatialDim);
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void test_polarity(const stk::mesh::EntityIdVector& nodeIds, const unsigned ordinal, const bool expectedPolarity)
  {
    stk::mesh::EntityVector nodes;
    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1u);
    EXPECT_TRUE(get_bulk().is_valid(elem1));

    for (auto nodeId : nodeIds)
    {
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
      EXPECT_TRUE(get_bulk().is_valid(node));
      nodes.push_back(node);
    }
    stk::mesh::EquivAndPositive result = is_side_equivalent_and_positive(get_bulk(), elem1, ordinal, nodes);
    EXPECT_TRUE(result.is_equiv);
    EXPECT_EQ(expectedPolarity, result.is_positive);
  }
};

typedef TestTextMesh TestQuad4;

TEST_F(TestQuad4, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  if (get_bulk().parallel_size() == 1)
  {
    inititalize_2D_mesh();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestQuad9;

TEST_F(TestQuad9, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,QUAD_9_2D,1,2,3,4,5,6,7,8,9";
  if (get_bulk().parallel_size() == 1)
  {
    inititalize_2D_mesh();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 5};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1, 5};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestTri3;

TEST_F(TestTri3, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,TRI_3_2D,1,2,3";
  if (get_bulk().parallel_size() == 1)
  {
    inititalize_2D_mesh();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestTri6;

TEST_F(TestTri6, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,TRI_6_2D,1,2,3,4,5,6";
  if (get_bulk().parallel_size() == 1)
  {
    inititalize_2D_mesh();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 4};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1, 4};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestHex8;

TEST_F(TestHex8, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 6, 5};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {5, 6, 2, 1};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestHex20;

TEST_F(TestHex20, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,HEX_20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 6, 5, 9, 14, 17, 13};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 5, 6, 2, 13, 17, 14, 9};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestPyramid5;

TEST_F(TestPyramid5, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 5};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {4, 3, 2, 1};
      unsigned ordinal = 4;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 5, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3, 4};
      unsigned ordinal = 4;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestPyramid13;

TEST_F(TestPyramid13, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,PYRAMID_13,1,2,3,4,5,6,7,8,9,10,11,12,13";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 5, 6, 11, 10};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 3, 2, 9, 8, 7, 6};
      unsigned ordinal = 4;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 5, 2, 10, 11, 6};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3, 4, 6, 7, 8, 9};
      unsigned ordinal = 4;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestTet4;

TEST_F(TestTet4, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,TET_4,1,2,3,4";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 4};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestTet10;

TEST_F(TestTet10, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,TET_10,1,2,3,4,5,6,7,8,9,10";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 4, 5, 9, 8};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 2, 8, 9, 5};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestWedge6;

TEST_F(TestWedge6, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,WEDGE_6,1,2,3,4,5,6";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 5, 4};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 3, 2};
      unsigned ordinal = 3;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 5, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3};
      unsigned ordinal = 3;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestWedge15;

TEST_F(TestWedge15, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,WEDGE_15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 5, 4, 7, 11, 13, 10};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 3, 2, 9, 8, 7};
      unsigned ordinal = 3;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 5, 2, 10, 13, 11, 7};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3, 7, 8, 9};
      unsigned ordinal = 3;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestShellLine2;

TEST_F(TestShellLine2, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,SHELL_LINE_2,1,2";
  if (get_bulk().parallel_size() == 1)
  {
    inititalize_2D_mesh();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestShellLine3;

TEST_F(TestShellLine3, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,SHELL_LINE_3,1,2,3";
  if (get_bulk().parallel_size() == 1)
  {
    inititalize_2D_mesh();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestShellTri3;

TEST_F(TestShellTri3, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,SHELL_TRI_3,1,2,3";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 3, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestShellTri6;

TEST_F(TestShellTri6, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,SHELL_TRI_6,1,2,3,4,5,6";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3, 4, 5, 6};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 3, 2, 6, 5, 4};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestShellQuad4;

TEST_F(TestShellQuad4, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,3,4";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3, 4};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 3, 2};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestShellQuad8;

TEST_F(TestShellQuad8, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,SHELL_QUAD_8,1,2,3,4,5,6,7,8";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3, 4, 5, 6, 7, 8};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 3, 2, 8, 7, 6, 5};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

typedef TestTextMesh TestShellQuad9;

TEST_F(TestShellQuad9, createNodeOrderingAndTestPolarity)
{
  std::string meshDesc = "0,1,SHELL_QUAD_9,1,2,3,4,5,6,7,8,9";
  if (get_bulk().parallel_size() == 1)
  {
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = true;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {2, 1, 3};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
    {
      stk::mesh::EntityIdVector nodeIds = {1, 4, 3, 2, 8, 7, 6, 5, 9};
      unsigned ordinal = 0;
      bool expectedPositivePolarity = false;
      test_polarity(nodeIds, ordinal, expectedPositivePolarity);
    }
  }
}

void run_one_parallel_boundary_polarity_test(stk::ParallelMachine comm,
                                             stk::mesh::BulkData::AutomaticAuraOption option)
{
  // Create bulk
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(3);
  builder.set_aura_option(option);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  stk::mesh::MetaData &meta = bulk->mesh_meta_data();

  // Register sideset updater
  if (!bulk->has_observer_type<stk::mesh::SidesetUpdater>()) {
    stk::mesh::Selector activeSelector = meta.universal_part();

    bulk->register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(*bulk, activeSelector));
  }

  // Create sideset parts
  std::vector<std::string> sidesetNames{"surface_1", "surface_block_1_QUAD4_1"};
  stk::mesh::PartVector sidesetParts;
  int sidesetId = 1;
  for(const std::string &name : sidesetNames)
  {
    stk::mesh::Part &part = meta.declare_part_with_topology(name, stk::topology::QUAD_4);
    meta.set_part_id(part, sidesetId);
    sidesetParts.push_back(&part);
  }

  meta.declare_part_subset(*sidesetParts[0], *sidesetParts[1]);

  // Build mesh
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
                         "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  bulk->initialize_face_adjacent_element_graph();
  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  // Create and populate sideset between elements 1 & 2, based on element 1
  stk::mesh::SideSet &sideset = bulk->create_sideset(*sidesetParts[0]);
  sideset.set_accept_all_internal_non_coincident_entries(false);

  std::vector<stk::mesh::EntityId> elementIds{1, 2};
  std::vector<stk::mesh::ConnectivityOrdinal> elementOrdinals{5, 4};
  std::vector<stk::mesh::Entity> elements(2);
  elements[0] = bulk->get_entity(stk::topology::ELEMENT_RANK, elementIds[0]);
  elements[1] = bulk->get_entity(stk::topology::ELEMENT_RANK, elementIds[1]);

  if (bulk->is_valid(elements[0])) {
    sideset.add({elements[0], elementOrdinals[0]});
  }

  // Setup surface to block mapping for the sideset updater
  stk::mesh::Part* block_1 = meta.get_part("block_1");
  EXPECT_TRUE(block_1 != nullptr);

  std::vector<const stk::mesh::Part*> touchingParts{block_1};
  meta.set_surface_to_block_mapping(sidesetParts[1], touchingParts);

  // Create faces
  bulk->modification_begin();
  int myProc = bulk->parallel_rank();
  stk::mesh::Entity face = bulk->declare_element_side(elements[myProc], elementOrdinals[myProc], sidesetParts);
  bulk->modification_end();

  // Make sure there is only 1 global face
  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(*bulk, counts);
  EXPECT_EQ(1u, counts[meta.side_rank()]);

  // Make sure it matches the one we created
  stk::mesh::EntityVector faces;
  stk::mesh::get_selected_entities(meta.universal_part(), bulk->buckets(meta.side_rank()), faces);
  EXPECT_EQ(1u, faces.size());
  EXPECT_EQ(face, faces[0]);


  // Get polarity results with free function code
  std::pair<bool,bool> freePolarity1 = stk::mesh::is_positive_sideset_polarity(*bulk, *sidesetParts[1], face,
                                                                               &meta.universal_part(), &sideset);
  std::pair<bool,bool> freePolarity2 = stk::mesh::is_positive_sideset_face_polarity(*bulk, face,
                                                                                    &meta.universal_part());

  // Get polarity results with util wrapper code
  stk::mesh::PolarityUtil helper(*bulk, meta.universal_part());
  std::pair<bool,bool> utilPolarity1 = helper.is_positive_sideset_polarity(*sidesetParts[1], face, &sideset);
  std::pair<bool,bool> utilPolarity2 = helper.is_positive_sideset_face_polarity(face);

  std::pair<bool,bool> expectedPolarity(true, true);

  EXPECT_EQ(expectedPolarity, freePolarity1);
  EXPECT_EQ(expectedPolarity, freePolarity2);

  EXPECT_EQ(expectedPolarity, utilPolarity1);
  EXPECT_EQ(expectedPolarity, utilPolarity2);
}

TEST(TestPolarity, oneParallelBoundaryBetweenTwoDifferentBlocks_withNoAura)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 2) {
    GTEST_SKIP();
  }

  run_one_parallel_boundary_polarity_test(comm, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(TestPolarity, oneParallelBoundaryBetweenTwoDifferentBlocks_withAura)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 2) {
    GTEST_SKIP();
  }

  run_one_parallel_boundary_polarity_test(comm, stk::mesh::BulkData::AUTO_AURA);
}

void run_two_parallel_boundaries_polarity_test(stk::ParallelMachine comm,
                                               stk::mesh::BulkData::AutomaticAuraOption option)
{
  // Create bulk
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(3);
  builder.set_aura_option(option);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  stk::mesh::MetaData &meta = bulk->mesh_meta_data();

  // Register sideset updater
  if (!bulk->has_observer_type<stk::mesh::SidesetUpdater>()) {
    stk::mesh::Selector activeSelector = meta.universal_part();

    bulk->register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(*bulk, activeSelector));
  }

  // Create sideset parts
  std::vector<std::string> sidesetNames{"surface_1", "surface_block_1_QUAD4_1"};
  stk::mesh::PartVector sidesetParts;
  int sidesetId = 1;
  for(const std::string &name : sidesetNames)
  {
    stk::mesh::Part &part = meta.declare_part_with_topology(name, stk::topology::QUAD_4);
    meta.set_part_id(part, sidesetId);
    sidesetParts.push_back(&part);
  }

  meta.declare_part_subset(*sidesetParts[0], *sidesetParts[1]);

  // Build mesh
  std::string meshDesc = "0,1,HEX_8, 1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                         "1,2,HEX_8, 5, 6, 7, 8, 9,10,11,12,block_2\n"
                         "0,3,HEX_8, 4, 3,14,13, 8, 7,16,15,block_1\n"
                         "0,4,HEX_8, 8, 7,16,15,12,11,18,17,block_1"
      ;
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,2, 1,0,2, 1,1,2, 0,1,2,
                                      0,2,0, 1,2,0, 0,2,1, 1,2,1,
                                      0,2,2, 1,2,2};

  bulk->initialize_face_adjacent_element_graph();
  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  // Create sideset
  stk::mesh::SideSet &sideset = bulk->create_sideset(*sidesetParts[0]);
  sideset.set_accept_all_internal_non_coincident_entries(false);

  // Populate sideset between elements 1 & 2, based on element 1
  std::vector<stk::mesh::EntityId> elementIds_face1{1, 2};
  std::vector<stk::mesh::ConnectivityOrdinal> elementOrdinals_face1{5, 4};
  std::vector<stk::mesh::Entity> elements_face1(2);
  elements_face1[0] = bulk->get_entity(stk::topology::ELEMENT_RANK, elementIds_face1[0]);
  elements_face1[1] = bulk->get_entity(stk::topology::ELEMENT_RANK, elementIds_face1[1]);

  if (bulk->is_valid(elements_face1[0])) {
    sideset.add({elements_face1[0], elementOrdinals_face1[0]});
  }

  // Populate sideset between elements 4 & 2, based on element 4
  std::vector<stk::mesh::EntityId> elementIds_face2{4, 2};
  std::vector<stk::mesh::ConnectivityOrdinal> elementOrdinals_face2{0, 2};
  std::vector<stk::mesh::Entity> elements_face2(2);
  elements_face2[0] = bulk->get_entity(stk::topology::ELEMENT_RANK, elementIds_face2[0]);
  elements_face2[1] = bulk->get_entity(stk::topology::ELEMENT_RANK, elementIds_face2[1]);

  if (bulk->is_valid(elements_face2[0])) {
    sideset.add({elements_face2[0], elementOrdinals_face2[0]});
  }

  // Setup surface to block mapping for the sideset updater
  stk::mesh::Part* block_1 = meta.get_part("block_1");
  EXPECT_TRUE(block_1 != nullptr);

  std::vector<const stk::mesh::Part*> touchingParts{block_1};
  meta.set_surface_to_block_mapping(sidesetParts[1], touchingParts);

  // Create faces
  bulk->modification_begin();
  int myProc = bulk->parallel_rank();
  stk::mesh::Entity face1 = bulk->declare_element_side(elements_face1[myProc], elementOrdinals_face1[myProc], sidesetParts);
  stk::mesh::Entity face2 = bulk->declare_element_side(elements_face2[myProc], elementOrdinals_face2[myProc], sidesetParts);
  bulk->modification_end();

  // Make sure there are only 2 global faces
  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(*bulk, counts);
  EXPECT_EQ(2u, counts[meta.side_rank()]);

  // Make sure it matches the one we created
  stk::mesh::EntityVector faces;
  stk::mesh::get_selected_entities(meta.universal_part(), bulk->buckets(meta.side_rank()), faces);
  EXPECT_EQ(2u, faces.size());
  EXPECT_EQ(face1, faces[0]);
  EXPECT_EQ(face2, faces[1]);

  {
    // Get polarity results for face1 with free function code
    std::pair<bool,bool> freePolarity1 = stk::mesh::is_positive_sideset_polarity(*bulk, *sidesetParts[1], face1,
                                                                                 &meta.universal_part(), &sideset);
    std::pair<bool,bool> freePolarity2 = stk::mesh::is_positive_sideset_face_polarity(*bulk, face1,
                                                                                      &meta.universal_part());

    // Get polarity results with util wrapper code
    stk::mesh::PolarityUtil helper(*bulk, meta.universal_part());
    std::pair<bool,bool> utilPolarity1 = helper.is_positive_sideset_polarity(*sidesetParts[1], face1, &sideset);
    std::pair<bool,bool> utilPolarity2 = helper.is_positive_sideset_face_polarity(face1);

    std::pair<bool,bool> expectedPolarity(true, true);

    EXPECT_EQ(expectedPolarity, freePolarity1);
    EXPECT_EQ(expectedPolarity, freePolarity2);

    EXPECT_EQ(expectedPolarity, utilPolarity1);
    EXPECT_EQ(expectedPolarity, utilPolarity2);
  }

  {
    // Get polarity results for face2 with free function code
    std::pair<bool,bool> freePolarity1 = stk::mesh::is_positive_sideset_polarity(*bulk, *sidesetParts[1], face2,
                                                                                 &meta.universal_part(), &sideset);
    std::pair<bool,bool> freePolarity2 = stk::mesh::is_positive_sideset_face_polarity(*bulk, face2,
                                                                                      &meta.universal_part());

    // Get polarity results with util wrapper code
    stk::mesh::PolarityUtil helper(*bulk, meta.universal_part());
    std::pair<bool,bool> utilPolarity1 = helper.is_positive_sideset_polarity(*sidesetParts[1], face2, &sideset);
    std::pair<bool,bool> utilPolarity2 = helper.is_positive_sideset_face_polarity(face2);

    std::pair<bool,bool> expectedPolarity(true, false);

    EXPECT_EQ(expectedPolarity, freePolarity1);
    EXPECT_EQ(expectedPolarity, freePolarity2);

    EXPECT_EQ(expectedPolarity, utilPolarity1);
    EXPECT_EQ(expectedPolarity, utilPolarity2);
  }
}

TEST(TestPolarity, twoParallelBoundariesBetweenTwoDifferentBlocks_withNoAura)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 2) {
    GTEST_SKIP();
  }

  run_two_parallel_boundaries_polarity_test(comm, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(TestPolarity, twoParallelBoundariesBetweenTwoDifferentBlocks_withAura)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 2) {
    GTEST_SKIP();
  }

  run_two_parallel_boundaries_polarity_test(comm, stk::mesh::BulkData::AUTO_AURA);
}
