#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>

namespace {

TEST(StkIo, sideset_with_bad_sideOrdinal)
{
  stk::ParallelMachine pm = stk::parallel_machine_world();
  const int p_size = stk::parallel_machine_size(pm);
  if(p_size != 1) { GTEST_SKIP(); }

  const std::string meshDesc =
      "textmesh:0,1,TET_4, 1,2,3,4, block_1 \
       |sideset:name=surface_1; data=1,1, 1,2, 1,3, 1,5";

  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(pm).create();
  EXPECT_ANY_THROW(stk::io::fill_mesh(meshDesc, *meshPtr));
}

TEST(StkIo, sideset_with_name_longer_than_32_chars)
{
  stk::ParallelMachine pm = stk::parallel_machine_world();
  const int p_size = stk::parallel_machine_size(pm);
  if(p_size != 1) { GTEST_SKIP(); }

  std::string partName("this_is_a_name_that_is_looooooooooooooooooonger_than_32_chars");
  const std::string meshDesc =
      "textmesh:0,1,TET_4, 1,2,3,4, block_1\n"
      "|sideset:name="+partName+"; data=1,1, 1,2, 1,3, 1,4";

  std::string exoFileName("output.exo");

  {
    std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(pm).create();
    stk::io::fill_mesh(meshDesc, *meshPtr);

    stk::io::write_mesh(exoFileName, *meshPtr);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(pm).create();
    stk::io::fill_mesh(exoFileName, *meshPtr);

    const stk::mesh::Part* part = meshPtr->mesh_meta_data().get_part(partName);
    ASSERT_TRUE(part != nullptr);
    EXPECT_LT(60u, part->name().size());
  }
}
}
