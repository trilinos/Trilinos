#include "gtest/gtest.h"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"
#include "stk_unit_test_utils/StkBalanceUnitTestSettings.hpp"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_balance/balance.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"

class BoundingBoxSearch : public stk::unit_test_util::MeshFixture
{
protected:
  void make_4_unit_quad_shell_connected_mesh()
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,7,6\n"
                 "0,2,SHELL_QUAD_4,2,3,8,7\n"
                 "0,3,SHELL_QUAD_4,3,4,9,8\n"
                 "0,4,SHELL_QUAD_4,4,5,10,9";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,7,6\n"
                 "0,2,SHELL_QUAD_4,2,3,8,7\n"
                 "1,3,SHELL_QUAD_4,3,4,9,8\n"
                 "1,4,SHELL_QUAD_4,4,5,10,9";
    }
    else if (numProcs == 3) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,7,6\n"
                 "0,2,SHELL_QUAD_4,2,3,8,7\n"
                 "1,3,SHELL_QUAD_4,3,4,9,8\n"
                 "2,4,SHELL_QUAD_4,4,5,10,9";
    }
    else if (numProcs == 4) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,7,6\n"
                 "1,2,SHELL_QUAD_4,2,3,8,7\n"
                 "2,3,SHELL_QUAD_4,3,4,9,8\n"
                 "3,4,SHELL_QUAD_4,4,5,10,9";
    }

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0,
      0,1,0, 1,1,0, 2,1,0, 3,1,0, 4,1,0
    };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void make_4_unit_quad_shell_single_gap_mesh(double gapSize)
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,8,7\n"
                 "0,2,SHELL_QUAD_4,2,3,9,8\n"
                 "0,3,SHELL_QUAD_4,4,5,11,10\n"
                 "0,4,SHELL_QUAD_4,5,6,12,11";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,8,7\n"
                 "0,2,SHELL_QUAD_4,2,3,9,8\n"
                 "1,3,SHELL_QUAD_4,4,5,11,10\n"
                 "1,4,SHELL_QUAD_4,5,6,12,11";
    }
    else if (numProcs == 3) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,8,7\n"
                 "0,2,SHELL_QUAD_4,2,3,9,8\n"
                 "1,3,SHELL_QUAD_4,4,5,11,10\n"
                 "2,4,SHELL_QUAD_4,5,6,12,11";
    }
    else if (numProcs == 4) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,8,7\n"
                 "1,2,SHELL_QUAD_4,2,3,9,8\n"
                 "2,3,SHELL_QUAD_4,4,5,11,10\n"
                 "3,4,SHELL_QUAD_4,5,6,12,11";
    }

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 2,0,0, 2+gapSize,0,0, 3+gapSize,0,0, 4+gapSize,0,0,
      0,1,0, 1,1,0, 2,1,0, 2+gapSize,1,0, 3+gapSize,1,0, 4+gapSize,1,0
    };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void make_4_unit_quad_shell_all_gap_mesh(double gapSize)
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,10,9\n"
                 "0,2,SHELL_QUAD_4,3,4,12,11\n"
                 "0,3,SHELL_QUAD_4,5,6,14,13\n"
                 "0,4,SHELL_QUAD_4,7,8,16,15";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,10,9\n"
                 "0,2,SHELL_QUAD_4,3,4,12,11\n"
                 "1,3,SHELL_QUAD_4,5,6,14,13\n"
                 "1,4,SHELL_QUAD_4,7,8,16,15";
    }
    else if (numProcs == 3) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,10,9\n"
                 "0,2,SHELL_QUAD_4,3,4,12,11\n"
                 "1,3,SHELL_QUAD_4,5,6,14,13\n"
                 "2,4,SHELL_QUAD_4,7,8,16,15";
    }
    else if (numProcs == 4) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,10,9\n"
                 "1,2,SHELL_QUAD_4,3,4,12,11\n"
                 "2,3,SHELL_QUAD_4,5,6,14,13\n"
                 "3,4,SHELL_QUAD_4,7,8,16,15";
    }

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 1+gapSize,0,0, 2+gapSize,0,0, 2+2*gapSize,0,0, 3+2*gapSize,0,0, 3+3*gapSize,0,0, 4+3*gapSize,0,0,
      0,1,0, 1,1,0, 1+gapSize,1,0, 2+gapSize,1,0, 2+2*gapSize,1,0, 3+2*gapSize,1,0, 3+3*gapSize,1,0, 4+3*gapSize,1,0
    };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void make_6_unit_tri_shell_connected_mesh()
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,5\n"
                 "0,2,SHELL_TRI_3,2,6,5\n"
                 "0,3,SHELL_TRI_3,2,3,6\n"
                 "0,4,SHELL_TRI_3,3,7,6\n"
                 "0,5,SHELL_TRI_3,3,4,7\n"
                 "0,6,SHELL_TRI_3,4,8,7";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,5\n"
                 "0,2,SHELL_TRI_3,2,6,5\n"
                 "0,3,SHELL_TRI_3,2,3,6\n"
                 "1,4,SHELL_TRI_3,3,7,6\n"
                 "1,5,SHELL_TRI_3,3,4,7\n"
                 "1,6,SHELL_TRI_3,4,8,7";
    }
    else if (numProcs == 3) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,5\n"
                 "0,2,SHELL_TRI_3,2,6,5\n"
                 "0,3,SHELL_TRI_3,2,3,6\n"
                 "1,4,SHELL_TRI_3,3,7,6\n"
                 "1,5,SHELL_TRI_3,3,4,7\n"
                 "2,6,SHELL_TRI_3,4,8,7";
    }
    else if (numProcs == 4) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,5\n"
                 "1,2,SHELL_TRI_3,2,6,5\n"
                 "1,3,SHELL_TRI_3,2,3,6\n"
                 "2,4,SHELL_TRI_3,3,7,6\n"
                 "2,5,SHELL_TRI_3,3,4,7\n"
                 "3,6,SHELL_TRI_3,4,8,7";
    }

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 2,0,0, 3,0,0,
      0,1,0, 1,1,0, 2,1,0, 3,1,0
    };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void make_6_unit_tri_shell_single_gap_mesh(double gapSize)
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,6\n"
                 "0,2,SHELL_TRI_3,2,7,6\n"
                 "0,3,SHELL_TRI_3,2,3,7\n"
                 "0,4,SHELL_TRI_3,4,9,8\n"
                 "0,5,SHELL_TRI_3,4,5,9\n"
                 "0,6,SHELL_TRI_3,5,10,9";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,6\n"
                 "0,2,SHELL_TRI_3,2,7,6\n"
                 "0,3,SHELL_TRI_3,2,3,7\n"
                 "1,4,SHELL_TRI_3,4,9,8\n"
                 "1,5,SHELL_TRI_3,4,5,9\n"
                 "1,6,SHELL_TRI_3,5,10,9";
    }
    else if (numProcs == 3) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,6\n"
                 "0,2,SHELL_TRI_3,2,7,6\n"
                 "0,3,SHELL_TRI_3,2,3,7\n"
                 "1,4,SHELL_TRI_3,4,9,8\n"
                 "1,5,SHELL_TRI_3,4,5,9\n"
                 "2,6,SHELL_TRI_3,5,10,9";
    }
    else if (numProcs == 4) {
      meshDesc = "0,1,SHELL_TRI_3,1,2,6\n"
                 "1,2,SHELL_TRI_3,2,7,6\n"
                 "1,3,SHELL_TRI_3,2,3,7\n"
                 "2,4,SHELL_TRI_3,4,9,8\n"
                 "2,5,SHELL_TRI_3,4,5,9\n"
                 "3,6,SHELL_TRI_3,5,10,9";
    }

    std::vector<double> coordinates;
    coordinates = {
      0,0,0, 1,0,0, 2,0,0,         2+gapSize,0,0, 3+gapSize,0,0,
      0,1,0, 1,1,0, 1+gapSize,1,0, 2+gapSize,1,0, 3+gapSize,1,0
    };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void make_4_unit_hex_single_gap_mesh(double gapSize)
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                 "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                 "0,3,HEX_8,13,14,15,16,17,18,19,20\n"
                 "0,4,HEX_8,17,18,19,20,21,22,23,24";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                 "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                 "1,3,HEX_8,13,14,15,16,17,18,19,20\n"
                 "1,4,HEX_8,17,18,19,20,21,22,23,24";
    }
    else if (numProcs == 3) {
      meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                 "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                 "1,3,HEX_8,13,14,15,16,17,18,19,20\n"
                 "2,4,HEX_8,17,18,19,20,21,22,23,24";
    }
    else if (numProcs == 4) {
      meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                 "1,2,HEX_8,5,6,7,8,9,10,11,12\n"
                 "2,3,HEX_8,13,14,15,16,17,18,19,20\n"
                 "3,4,HEX_8,17,18,19,20,21,22,23,24";
    }

    std::vector<double> coordinates = {
      0,0,0,         0,1,0,         0,1,1,         0,0,1,
      1,0,0,         1,1,0,         1,1,1,         1,0,1,
      2,0,0,         2,1,0,         2,1,1,         2,0,1,
      2+gapSize,0,0, 2+gapSize,1,0, 2+gapSize,1,1, 2+gapSize,0,1,
      3+gapSize,0,0, 3+gapSize,1,0, 3+gapSize,1,1, 3+gapSize,0,1,
      4+gapSize,0,0, 4+gapSize,1,0, 4+gapSize,1,1, 4+gapSize,0,1
    };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void make_two_particle_mesh()
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,PARTICLE,1\n"
                 "0,2,PARTICLE,2";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,PARTICLE,1\n"
                 "1,2,PARTICLE,2";
    }
    std::vector<double> coordinates = { 0,0,0, 2,0,0 };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void make_particle_unit_quad_shell_mesh()
  {
    const int numProcs = get_parallel_size();

    std::string meshDesc;
    if (numProcs == 1) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,3,4\n"
                 "0,2,PARTICLE,5";
    }
    else if (numProcs == 2) {
      meshDesc = "0,1,SHELL_QUAD_4,1,2,3,4\n"
                 "1,2,PARTICLE,5";
    }
    std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 2,0,0 };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  bool check_edges(const std::vector<stk::balance::GraphEdge> & graphEdges,
                   const std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId>> & expectedEdges)
  {
    stk::mesh::BulkData & bulk = get_bulk();
    bool match = true;

    const size_t globalNumEdges = stk::get_global_sum(bulk.parallel(), graphEdges.size());
    if ((bulk.parallel_rank() == 0) && (globalNumEdges != expectedEdges.size())) {
      std::cout << "Actual size (" << globalNumEdges << ") not equal to expected size ("
                << expectedEdges.size() << ")" << std::endl;
      match = false;
    }

    for (const auto & expected : expectedEdges) {
      stk::mesh::Entity firstEntity = bulk.get_entity(stk::topology::ELEM_RANK, expected.first);
      if (!(bulk.is_valid(firstEntity) && bulk.bucket(firstEntity).owned())) continue;

      const auto foundEdge = std::find_if(graphEdges.begin(), graphEdges.end(),
          [&](const stk::balance::GraphEdge & curr){
              return ( (bulk.identifier(curr.vertex1()) == expected.first) &&
                       (curr.vertex2_id() == expected.second) );
          });

      if (foundEdge == graphEdges.end()) {
        std::cout << "Missing graph edge: (" << expected.first << ", "
                  << expected.second << ")" << std::endl;
        match = false;
      }
    }

    return stk::is_true_on_all_procs(bulk.parallel(), match);
  }

};



TEST_F(BoundingBoxSearch, fourUnitQuadShell_connected)
{
//       6-----7-----8-----9----10
//       |     |     |     |     |
//       |  1  |  2  |  3  |  4  |
//       |     |     |     |     |
//       1-----2-----3-----4-----5
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.7);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_4_unit_quad_shell_connected_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {}));
}

TEST_F(BoundingBoxSearch, fourUnitQuadShell_tooLargeTolerance)
{
//       6-----7-----8-----9----10
//       |     |     |     |     |
//       |  1  |  2  |  3  |  4  |
//       |     |     |     |     |
//       1-----2-----3-----4-----5
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(1.1);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_4_unit_quad_shell_connected_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{1,3}, {3,1}, {2,4}, {4,2}}));
}

TEST_F(BoundingBoxSearch, fourUnitQuadShell_smallGap)
{
//       7-----8-----9 10----11---12
//       |     |     | |     |     |
//       |  1  |  2  | |  3  |  4  |
//       |     |     | |     |     |
//       1-----2-----3 4-----5-----6
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.7);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  const double gapSize = 0.1;
  make_4_unit_quad_shell_single_gap_mesh(gapSize);

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{2,3}, {3,2}}));
}

TEST_F(BoundingBoxSearch, fourUnitQuadShell_largeGap)
{
//       7-----8-----9     10----11---12
//       |     |     |     |     |     |
//       |  1  |  2  |     |  3  |  4  |
//       |     |     |     |     |     |
//       1-----2-----3     4-----5-----6
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.6);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  const double gapSize = 1.0;
  make_4_unit_quad_shell_single_gap_mesh(gapSize);

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {}));
}

TEST_F(BoundingBoxSearch, fourUnitQuadShell_allGaps)
{
//       9----10 11---12 13---14 15---16
//       |     | |     | |     | |     |
//       |  1  | |  2  | |  3  | |  4  |
//       |     | |     | |     | |     |
//       1-----2 3-----4 5-----6 7-----8
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.7);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  const double gapSize = 0.1;
  make_4_unit_quad_shell_all_gap_mesh(gapSize);

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{1,2}, {2,1}, {2,3}, {3,2}, {3,4}, {4,3}}));
}

TEST_F(BoundingBoxSearch, sixUnitTriShell_connected)
{
//       5-----6-----7-----8
//       |\  2 |\  4 |\  6 |
//       |  \  |  \  |  \  |
//       | 1  \| 3  \| 5  \|
//       1-----2-----3-----4
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.7);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_6_unit_tri_shell_connected_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {}));
}

TEST_F(BoundingBoxSearch, sixUnitTriShell_smallGap)
{
//       6-----7 8-----9----10
//       |\  2 |\ \  4 |\  6 |
//       |  \  |  \ \  |  \  |
//       | 1  \| 3  \ \| 5  \|
//       1-----2-----3 4-----5
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.09);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  const double gapSize = 0.1;
  make_6_unit_tri_shell_single_gap_mesh(gapSize);

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{3,4}, {4,3}}));
}

TEST_F(BoundingBoxSearch, sixUnitTriShell_largeGap)
{
//       6-----7     8-----9----10
//       |\  2 |\     \  4 |\  6 |
//       |  \  |  \     \  |  \  |
//       | 1  \| 3  \     \| 5  \|
//       1-----2-----3     4-----5
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.1);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  const double gapSize = 1.0;
  make_6_unit_tri_shell_single_gap_mesh(gapSize);

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {}));
}

TEST_F(BoundingBoxSearch, fourUnitHex_smallGap)
{
//           3--------------7-------------11      15-------------19------------23
//          /|             /|             /|     /|             /|             /|
//         / |            / |            / |    / |            / |            / |
//        /  |           /  |           /  |   /  |           /  |           /  |
//       4--------------8-------------12   |  16-------------20------------24   |
//       |   |          |   |          |   |  |   |          |   |          |   |
//       |   |     1    |   |     2    |   |  |   |     3    |   |     4    |   |
//       |   |          |   |          |   |  |   |          |   |          |   |
//       |   2----------|---6----------|--10  |   14---------|---18---------|--22
//       |  /           |  /           |  /   |  /           |  /           |  /
//       | /            | /            | /    | /            | /            | /
//       |/             |/             |/     |/             |/             |/
//       1--------------5--------------9      13-------------17------------21
//  z
//  | y
//  |/
//  o---x
//
  if (stk::parallel_machine_size(get_comm()) > 4) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.7);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  const double gapSize = 0.1;
  make_4_unit_hex_single_gap_mesh(gapSize);

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{2,3}, {3,2}}));
}

TEST_F(BoundingBoxSearch, particleParticle_smallTolerance)
{
//
//      o             o
//       (1)           (2)
//
// (Particle spacing: 2 units)
//
  if (stk::parallel_machine_size(get_comm()) > 2) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForParticleSearch(2.1);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_two_particle_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {}));
}

TEST_F(BoundingBoxSearch, particleParticle_largeTolerance)
{
//
//      o             o
//       (1)           (2)
//
// (Particle spacing: 2 units)
//
  if (stk::parallel_machine_size(get_comm()) > 2) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForParticleSearch(4.1);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_two_particle_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{1,2}, {2,1}}));
}

TEST_F(BoundingBoxSearch, particleQuadShell_smallTolerance)
{
//       4-----3
//       |     |
//       |  1  |
//       |     |
//       1-----2     o (Node 5, Elem 2)
//
  if (stk::parallel_machine_size(get_comm()) > 2) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.6);
  balanceSettings.setToleranceForParticleSearch(1.2);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_particle_unit_quad_shell_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {}));
}

TEST_F(BoundingBoxSearch, particleQuadShell_smallParticleTolerance)
{
//       4-----3
//       |     |
//       |  1  |
//       |     |
//       1-----2     o (Node 5, Elem 2)
//
  if (stk::parallel_machine_size(get_comm()) > 2) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(1.1);
  balanceSettings.setToleranceForParticleSearch(1.2);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_particle_unit_quad_shell_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{1,2}, {2,1}}));
}

TEST_F(BoundingBoxSearch, particleQuadShell_smallFaceTolerance)
{
//       4-----3
//       |     |
//       |  1  |
//       |     |
//       1-----2     o (Node 5, Elem 2)
//
  if (stk::parallel_machine_size(get_comm()) > 2) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(0.6);
  balanceSettings.setToleranceForParticleSearch(2.2);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_particle_unit_quad_shell_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{1,2}, {2,1}}));
}

TEST_F(BoundingBoxSearch, particleQuadShell_largeTolerance)
{
//       4-----3
//       |     |
//       |  1  |
//       |     |
//       1-----2     o (Node 5, Elem 2)
//
  if (stk::parallel_machine_size(get_comm()) > 2) return;

  stk::unit_test_util::StkBalanceUnitTestSettings balanceSettings;
  balanceSettings.setToleranceForFaceSearch(1.1);
  balanceSettings.setToleranceForParticleSearch(2.2);
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());

  make_particle_unit_quad_shell_mesh();

  std::vector<stk::balance::GraphEdge> graphEdges;
  stk::balance::internal::addGraphEdgesUsingBBSearch(get_bulk(),
                                                     balanceSettings,
                                                     graphEdges,
                                                     get_meta().universal_part());

  EXPECT_TRUE(check_edges(graphEdges, {{1,2}, {2,1}}));
}
