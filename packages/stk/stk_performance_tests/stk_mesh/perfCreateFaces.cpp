#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>

namespace
{

size_t get_num_global_faces(const stk::mesh::BulkData &bulk)
{
  std::vector<size_t> meshCounts;
  stk::mesh::comm_mesh_counts(bulk, meshCounts);
  return meshCounts[stk::topology::FACE_RANK];
}

//==============================================================================
class CreateFacesClassicPerformance : public stk::unit_test_util::PerformanceTester
{
public:
  CreateFacesClassicPerformance(stk::mesh::BulkData &bulk)
    : stk::unit_test_util::PerformanceTester(bulk.parallel()),
      bulkData(bulk),
      locallyOwnedElements(bulk.mesh_meta_data().locally_owned_part())
  {
  }

protected:
  virtual void run_algorithm_to_time()
  {
    stk::mesh::create_faces(bulkData, locallyOwnedElements);
  }

  virtual size_t get_value_to_output_as_iteration_count()
  {
    return get_num_global_faces(bulkData);
  }

  stk::mesh::BulkData &bulkData;
  stk::mesh::Selector locallyOwnedElements;
};

//==============================================================================
class CreateFacesPerformance : public stk::unit_test_util::PerformanceTester
{
public:
  CreateFacesPerformance(stk::mesh::BulkData &bulk) :
    stk::unit_test_util::PerformanceTester(bulk.parallel()),
    bulkData(bulk),
    locallyOwnedElements(bulk.mesh_meta_data().locally_owned_part())
  {
  }

protected:
  virtual void run_algorithm_to_time()
  {
    stk::mesh::experimental::create_faces(bulkData, locallyOwnedElements);
  }

  virtual size_t get_value_to_output_as_iteration_count()
  {
    return get_num_global_faces(bulkData);
  }

  stk::mesh::BulkData &bulkData;
  stk::mesh::Selector locallyOwnedElements;
};

//==============================================================================
class CreateFacesClassicPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
  void run_create_faces_perf_test()
  {
    CreateFacesClassicPerformance perfTester(get_bulk());
    perfTester.run_performance_test();
  }

  std::string get_mesh_spec()
  {
    return stk::unit_test_util::get_option("-file", "NO_FILE_SPECIFIED");
  }
};

//==============================================================================
class CreateFacesPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
  void run_create_faces_perf_test()
  {
    CreateFacesPerformance perfTester(get_bulk());
    perfTester.run_performance_test();
  }

  std::string get_mesh_spec()
  {
    return stk::unit_test_util::get_option("-file", "NO_FILE_SPECIFIED");
  }
};

//==============================================================================
TEST_F(CreateFacesClassicPerformanceTest, read_mesh)
{
  setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

  run_create_faces_perf_test();
}

TEST_F(CreateFacesPerformanceTest, read_mesh_with_auto_decomp)
{
  allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
  stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

  run_create_faces_perf_test();
}

TEST_F(CreateFacesPerformanceTest, read_mesh)
{
  setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

  run_create_faces_perf_test();
}

}

