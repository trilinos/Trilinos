#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_io/IossBridge.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SideSetUtil.hpp>
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/IOHelpers.hpp>
#include "stk_unit_test_utils/FaceTestingUtils.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "IOMeshFixture.hpp"
#include <stk_unit_test_utils/BuildMesh.hpp>
#include "Ioss_SideSet.h"
#include "Ioss_SideBlock.h"

namespace {
using stk::unit_test_util::build_mesh_no_simple_fields;

struct SideSetDescription
{
  SideSetDescription(bool should_use_all_face_sides_ = false, int side_block_dim_=0, int side_set_dim_=0, stk::mesh::EntityRank side_block_rank_=stk::topology::END_RANK) :
    should_use_all_face_sides(should_use_all_face_sides_),
    side_block_dim(side_block_dim_),
    side_set_dim(side_set_dim_),
    side_block_rank(side_block_rank_)
  {}

  bool should_use_all_face_sides;
  int side_block_dim;
  int side_set_dim;
  stk::mesh::EntityRank side_block_rank;
};

bool operator==(const SideSetDescription& lhs, const SideSetDescription& rhs)
{
  return lhs.should_use_all_face_sides == rhs.should_use_all_face_sides &&
         lhs.side_block_dim == rhs.side_block_dim &&
         lhs.side_set_dim   == rhs.side_set_dim &&
         lhs.side_block_rank == rhs.side_block_rank;
}

std::ostream& operator<<(std::ostream& os, const SideSetDescription& desc)
{
  os << "should use all_face_sides = " << desc.should_use_all_face_sides 
     << ",  block dim = " << desc.side_block_dim 
     << ", set dim = " << desc.side_set_dim
     << ", block rank = " << desc.side_block_rank;

  return os;
}

SideSetDescription get_sideset_desc(const std::string& meshDesc, unsigned spatialDim, bool enableAllFaceSides)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::io::StkMeshIoBroker broker(bulk->parallel());
  broker.set_enable_all_face_sides_shell_topo(enableAllFaceSides);
  stk::io::fill_mesh_preexisting(broker, meshDesc, *bulk);

  std::shared_ptr<Ioss::Region> region = broker.get_input_ioss_region();
  const Ioss::SideSetContainer& side_sets = region->get_sidesets();
  Ioss::SideSet* side_set = *(side_sets.begin());
  Ioss::SideBlock* side_block =  *(side_set->get_side_blocks().begin());

  return SideSetDescription(stk::io::should_use_all_face_sides(side_block),
                            stk::io::get_max_par_dimension(side_block),
                            stk::io::get_max_par_dimension(side_set),
                            stk::io::get_side_rank(side_block));
}


void check_should_use_all_face_sides(const std::string& meshDesc, unsigned spatialDim, bool enableAllFaceSides,
                                     bool val_for_element_block, bool val_for_side_block)
{

  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::io::StkMeshIoBroker broker(bulk->parallel());
  broker.set_enable_all_face_sides_shell_topo(enableAllFaceSides);
  stk::io::fill_mesh_preexisting(broker, meshDesc, *bulk);

  std::shared_ptr<Ioss::Region> region = broker.get_input_ioss_region();
  const Ioss::ElementBlockContainer& element_blocks = region->get_element_blocks();
  Ioss::ElementBlock* block = *(element_blocks.begin());
  const Ioss::SideSetContainer& side_sets = region->get_sidesets();
  Ioss::SideSet* side_set = *(side_sets.begin());
  Ioss::SideBlock* side_block =  *(side_set->get_side_blocks().begin());

  EXPECT_EQ(stk::io::should_use_all_face_sides(block), val_for_element_block);
  EXPECT_EQ(stk::io::should_use_all_face_sides(side_block), val_for_side_block);
  EXPECT_EQ(stk::io::should_use_all_face_sides(static_cast<Ioss::EntityBlock*>(side_block)), val_for_side_block);
}

}


TEST(SidesetDimension, Shells)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "textmesh: 0,1,SHELL_QUAD_4,1,2,3,4 |coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0 |sideset:name=surface_1; data=1,3; split=topology";
  bool enableAllFaceSides = false;
  unsigned spatialDim = 3;

  check_should_use_all_face_sides(meshDesc, spatialDim, enableAllFaceSides, false, false);
  EXPECT_EQ(get_sideset_desc(meshDesc, spatialDim, enableAllFaceSides), SideSetDescription(false, 1, 1, stk::topology::EDGE_RANK));

  enableAllFaceSides = true;
  check_should_use_all_face_sides(meshDesc, spatialDim, enableAllFaceSides, true, true);
  EXPECT_EQ(get_sideset_desc(meshDesc, spatialDim, enableAllFaceSides), SideSetDescription(true, 2, 2, stk::topology::FACE_RANK));
}

TEST(SidesetDimension, TwoD)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "textmesh: 0,1,QUAD_4_2D,1,2,3,4 |dimension: 2 |coordinates: 0,0, 1,0, 1,1, 0,1, |sideset:name=surface_1; data=1,3; split=topology";
  bool enableAllFaceSides = false;
  unsigned spatialDim = 2;

  check_should_use_all_face_sides(meshDesc, spatialDim, enableAllFaceSides, false, false);
  EXPECT_EQ(get_sideset_desc(meshDesc, spatialDim, enableAllFaceSides), SideSetDescription(false, 1, 1, stk::topology::EDGE_RANK));

  enableAllFaceSides = true;
  check_should_use_all_face_sides(meshDesc, spatialDim, enableAllFaceSides, false, false);
  EXPECT_EQ(get_sideset_desc(meshDesc, spatialDim, enableAllFaceSides), SideSetDescription(false, 1, 1, stk::topology::EDGE_RANK));

}

TEST(SidesetDimension, Hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "generated:2x2x2|sideset:x";
  bool enableAllFaceSides = false;
  unsigned spatialDim = 3;

  check_should_use_all_face_sides(meshDesc, spatialDim, enableAllFaceSides, false, false);
  EXPECT_EQ(get_sideset_desc(meshDesc, spatialDim, enableAllFaceSides), SideSetDescription(false, 2, 2, stk::topology::FACE_RANK));

  enableAllFaceSides = true;
  check_should_use_all_face_sides(meshDesc, spatialDim, enableAllFaceSides, false, false);
  EXPECT_EQ(get_sideset_desc(meshDesc, spatialDim, enableAllFaceSides), SideSetDescription(false, 2, 2, stk::topology::FACE_RANK));

}

TEST(SidesetRank, ShellsWithAllFaceSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "textmesh: 0,1,SHELL_QUAD_4,1,2,3,4 |coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0 |sideset:name=surface_1; data=1,3; split=topology";
  unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::io::StkMeshIoBroker broker(bulk->parallel());
  broker.set_enable_all_face_sides_shell_topo(true);
  stk::io::fill_mesh_preexisting(broker, meshDesc, *bulk);

  EXPECT_EQ(bulk->mesh_meta_data().get_part("SURFACE_1")->primary_entity_rank(), stk::topology::FACE_RANK);
}

TEST(SidesetRank, ShellsWithoutAllFaceSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "textmesh: 0,1,SHELL_QUAD_4,1,2,3,4 |coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0 |sideset:name=surface_1; data=1,3; split=topology";
  unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::io::StkMeshIoBroker broker(bulk->parallel());
  broker.set_enable_all_face_sides_shell_topo(false);
  stk::io::fill_mesh_preexisting(broker, meshDesc, *bulk);

  EXPECT_EQ(bulk->mesh_meta_data().get_part("SURFACE_1")->primary_entity_rank(), stk::topology::EDGE_RANK);
}