#include "gtest/gtest.h"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_io/FillMesh.hpp>

namespace
{

TEST(SelfContainedModCycle, create_all_sides_throw)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::io::fill_mesh("generated:1x1x2", bulk);

  bulk.modification_begin(); 
  EXPECT_ANY_THROW(stk::mesh::create_all_sides(bulk, meta.universal_part(), stk::mesh::PartVector{}, true));
}

TEST(SelfContainedModCycle, create_interior_block_boundary_sides_throw)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::io::fill_mesh("generated:1x1x2", bulk);

  bulk.modification_begin(); 
  EXPECT_ANY_THROW(stk::mesh::create_interior_block_boundary_sides(bulk, meta.universal_part(),
                                                                   stk::mesh::PartVector{}));
}

TEST(SelfContainedModCycle, create_exposed_block_boundary_sides_throw)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::io::fill_mesh("generated:1x1x2", bulk);

  bulk.modification_begin(); 
  EXPECT_ANY_THROW(stk::mesh::create_exposed_block_boundary_sides(bulk, meta.universal_part(),
                                                                  stk::mesh::PartVector{}));
}

TEST(SelfContainedModCycle, create_all_block_boundary_sides_throw)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::io::fill_mesh("generated:1x1x2", bulk);

  bulk.modification_begin(); 
  EXPECT_ANY_THROW(stk::mesh::create_all_block_boundary_sides(bulk, meta.universal_part(),
                                                              stk::mesh::PartVector{}));
}

TEST(SelfContainedModCycle, create_exposed_block_boundary_sides_selector_throw)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::io::fill_mesh("generated:1x1x2", bulk);

  bulk.modification_begin(); 
  EXPECT_ANY_THROW(stk::mesh::create_exposed_block_boundary_sides(bulk, meta.universal_part(),
                                                                  stk::mesh::PartVector{}, meta.universal_part()));
}

}
