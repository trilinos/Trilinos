// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>
#include <stddef.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <string>
#include <vector>
#include "mpi.h"
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
namespace stk { namespace mesh { struct MeshIndex; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace
{
using stk::unit_test_util::build_mesh;

TEST(ForEachEntityRunNoThreads, test_counting_nodes)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) != 2) { GTEST_SKIP(); }

  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkData = *bulkPtr;

  std::string generatedMeshSpec = "generated:1x1x4";
  stk::io::fill_mesh(generatedMeshSpec, bulkData);

  {//version that uses selector and functor that takes MeshIndex
    unsigned numNodes = 0;
    stk::mesh::for_each_entity_run_no_threads(bulkData, stk::topology::NODE_RANK, metaData.universal_part(),
      [&numNodes](const stk::mesh::BulkData & mesh, const stk::mesh::MeshIndex & meshIndex)
      {
        stk::mesh::Entity entity = stk::mesh::impl::get_entity(meshIndex);
        if(mesh.is_valid(entity)) {
          numNodes++;
        }
      }
    );
    EXPECT_EQ(16u, numNodes);
  }
  {//version that uses selector and functor that takes Entity
    unsigned numNodes = 0;
    stk::mesh::for_each_entity_run_no_threads(bulkData, stk::topology::NODE_RANK, metaData.universal_part(),
      [&numNodes](const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity)
      {
        if(mesh.is_valid(entity)) {
          numNodes++;
        }
      }
    );
    EXPECT_EQ(16u, numNodes);
  }
  {//version with no selector and functor that takes MeshIndex
    unsigned numNodes = 0;
    stk::mesh::for_each_entity_run_no_threads(bulkData, stk::topology::NODE_RANK,
      [&numNodes](const stk::mesh::BulkData & mesh, const stk::mesh::MeshIndex & meshIndex)
      {
        stk::mesh::Entity entity = stk::mesh::impl::get_entity(meshIndex);
        if(mesh.is_valid(entity)) {
          numNodes++;
        }
      }
    );
    EXPECT_EQ(16u, numNodes);
  }
  {//version with no selector and functor that takes Entity
    unsigned numNodes = 0;
    stk::mesh::for_each_entity_run_no_threads(bulkData, stk::topology::NODE_RANK,
      [&numNodes](const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity)
      {
        if(mesh.is_valid(entity)) {
          numNodes++;
        }
      }
    );
    EXPECT_EQ(16u, numNodes);
  }
}

TEST(SettingLocalIds, bucket_loop_iterates_in_same_order_as_for_each_entity_run)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) != 2) { GTEST_SKIP(); }

  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkData = *bulkPtr;

  std::string generatedMeshSpec = "generated:1x1x4";
  stk::io::fill_mesh(generatedMeshSpec, bulkData);

  const unsigned numNodes = stk::mesh::count_selected_entities(metaData.universal_part(),
                                                               bulkData.buckets(stk::topology::NODE_RANK));
  ASSERT_EQ(16u, numNodes);
  std::vector<unsigned> localIds(numNodes);
  unsigned index = 0;
  const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, metaData.universal_part());
  for(const stk::mesh::Bucket* bptr : buckets) {
    for(stk::mesh::Entity node : *bptr) {
      if (bulkData.is_valid(node)) {
        bulkData.set_local_id(node, index);
        ++index;
      }
    }
  }
  stk::mesh::for_each_entity_run(bulkData, stk::topology::NODE_RANK,
    [&localIds](const stk::mesh::BulkData & mesh, const stk::mesh::MeshIndex & meshIndex)
    {
      stk::mesh::Entity node = stk::mesh::impl::get_entity(meshIndex);
      if(mesh.is_valid(node))
      {
        unsigned localId = mesh.local_id(node);
        localIds[localId] = localId;
      }
    }
  );
  for(size_t i=0; i<numNodes; i++)
  {
    EXPECT_EQ(i, localIds[i]);
  }
}

}
