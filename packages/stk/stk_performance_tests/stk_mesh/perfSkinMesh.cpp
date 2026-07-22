#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/DestroyRelations.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
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

class SkinMeshPerformance : public stk::unit_test_util::PerformanceTester
{
public:
  SkinMeshPerformance(stk::mesh::BulkData &bulk)
    : stk::unit_test_util::PerformanceTester(bulk.parallel()),
      bulkData(bulk),
      thingToSkin(bulk.mesh_meta_data().universal_part()),
      skinPart(bulk.mesh_meta_data().declare_part("skinPart"))
  {
  }

  void delete_sides()
  {
    const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, meta.locally_owned_part(), elems);
    bulkData.modification_begin();

    for(stk::mesh::Entity elem : elems)
    {
      stk::mesh::destroy_relations(bulkData, elem, meta.side_rank());
    }

    stk::mesh::EntityVector sides;
    stk::mesh::get_entities(bulkData, meta.side_rank(), meta.locally_owned_part(), sides);
    for(stk::mesh::Entity side : sides)
    {
      bulkData.destroy_entity(side);
    }

    bulkData.modification_end();
  }


protected:
  virtual void run_algorithm_to_time()
  {
    stk::mesh::skin_mesh(bulkData, thingToSkin, {&skinPart});
  }

  virtual size_t get_value_to_output_as_iteration_count()
  {
    return get_num_global_faces(bulkData);
  }

  stk::mesh::BulkData &bulkData;
  stk::mesh::Selector thingToSkin;
  stk::mesh::Part &skinPart;
};

class SkinMeshPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
  void run_skin_mesh_perf_test()
  {
    SkinMeshPerformance perfTester(get_bulk());
    const int numTimes = 20;
    for (int i = 0; i < numTimes; ++i) {
      perfTester.delete_sides();
      perfTester.run_performance_test();
    }
  }
  std::string get_mesh_spec()
  {
    return stk::unit_test_util::get_option("-file", "NO_FILE_SPECIFIED");
  }
};

TEST_F(SkinMeshPerformanceTest, read_mesh_with_auto_decomp)
{
  allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
  stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

  run_skin_mesh_perf_test();
}

TEST_F(SkinMeshPerformanceTest, read_mesh)
{
  setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

  run_skin_mesh_perf_test();
}

}

