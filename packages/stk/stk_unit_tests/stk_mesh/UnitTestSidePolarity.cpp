#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>


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
    allocate_meta(spatialDim);
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
