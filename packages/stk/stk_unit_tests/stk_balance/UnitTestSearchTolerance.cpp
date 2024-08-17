#include <stk_mesh/base/Field.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/setup/DefaultSettings.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/search_tolerance_algs/SecondShortestEdgeFaceSearchTolerance.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>

namespace
{

class SearchToleranceTest : public stk::unit_test_util::MeshFixture {};

TEST_F(SearchToleranceTest, faceOfCube)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  if (numProcs == 1) {
    setup_mesh("generated:1x1x1|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);

    const stk::mesh::FieldBase* coordField = get_bulk().mesh_meta_data().coordinate_field();
    stk::balance::SecondShortestEdgeFaceSearchTolerance faceSearchTolerance;

    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::Entity face = get_bulk().begin_faces(elem)[0];

    double searchTol = faceSearchTolerance.compute(get_bulk(), *coordField, get_bulk().begin_nodes(face), get_bulk().num_nodes(face));

    const double epsilon = 1.e-6;
    EXPECT_NEAR(0.15, searchTol, epsilon);
  }
}

void change_position_of_node_1(stk::mesh::BulkData& bulk)
{
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  const stk::mesh::FieldBase* coordField = bulk.mesh_meta_data().coordinate_field();
  double* coordData = static_cast<double*>(stk::mesh::field_data(*coordField, node1));
  coordData[0] = 0.25;
  coordData[1] = 0.0;
  coordData[2] = 0.5;
}

TEST_F(SearchToleranceTest, faceWithDifferentEdgeLengths)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  if (numProcs == 1) {
    setup_mesh("generated:1x1x1|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);

    change_position_of_node_1(get_bulk());

    stk::balance::SecondShortestEdgeFaceSearchTolerance faceSearchTolerance;

    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::Entity face = get_bulk().begin_faces(elem)[0];
    const stk::mesh::FieldBase* coordField = get_bulk().mesh_meta_data().coordinate_field();

    double searchTol = faceSearchTolerance.compute(get_bulk(), *coordField, get_bulk().begin_nodes(face), get_bulk().num_nodes(face));

    double secondShortestEdgeLength = 0.9013878;
    const double epsilon = 1.e-6;
    EXPECT_NEAR(0.15*secondShortestEdgeLength, searchTol, epsilon);
  }
}

class SearchToleranceTester : public stk::unit_test_util::MeshFixture
{
protected:

  void make_two_separated_hex_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16";
    double eps = 0.1;
    std::vector<double> coordinates {
      0,0,0,
      1,0,0,
      1,0,1,
      0,0,1,
      0,1,0,
      1,1,0,
      1,1,1,
      0,1,1,

      0,eps+1,0,
          1,eps+1,0,
          1,eps+1,1,
          0,eps+1,1,
          0,eps+2,0,
          1,eps+2,0,
          1,eps+2,1,
          0,eps+2,1,
    };
    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  unsigned get_num_search_results_with_app_settings(const stk::balance::GraphCreationSettings &balanceSettings)
  {
    make_two_separated_hex_mesh();
    stk::mesh::Selector thingsToSearch = get_meta().universal_part();
    stk::balance::internal::SearchElemPairs searchResults = stk::balance::internal::getBBIntersectionsForFacesParticles(get_bulk(), balanceSettings, thingsToSearch);
    return searchResults.size();
  }
};

TEST_F(SearchToleranceTester, constantTolerance)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    stk::balance::GraphCreationSettings balanceSettings;
    balanceSettings.setToleranceForFaceSearch(stk::balance::DefaultSettings::faceSearchAbsTol);
    const unsigned numSelfInteractions = 2;
    EXPECT_EQ(numSelfInteractions, get_num_search_results_with_app_settings(balanceSettings));
  }
}

TEST_F(SearchToleranceTester, secondShortestEdgeFaceSearchTolerance)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    stk::balance::GraphCreationSettings balanceSettings;
    const unsigned numSelfPlusSymmetricInteractions = 4;
    EXPECT_EQ(numSelfPlusSymmetricInteractions, get_num_search_results_with_app_settings(balanceSettings));
  }
}

}
