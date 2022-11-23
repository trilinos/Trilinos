#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include "stk_util/parallel/Parallel.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"

namespace
{

TEST(StkIo, readingParallelFilesMissingParallelCommInfo)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(communicator).create();
        // file from Salinas output is missing parallel info
        EXPECT_THROW(stk::io::fill_mesh("twoHexMissingParallelInfo.e", *bulk), std::runtime_error);
    }
}


TEST(StkIo, checkCanonicalName)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) != 1) return;

  stk::io::StkMeshIoBroker stkIo;
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(communicator).create();

  std::string meshSpec = "textmesh:0,1,QUAD_4_2D,1,2,3,4,Unspecified-2-QUAD|dimension:2";
  std::string partName = "Unspecified-2-QUAD";

  stk::io::fill_mesh(meshSpec, *bulk, stkIo);

  stk::mesh::Part* metaPart = bulk->mesh_meta_data().get_part(partName);
  stk::mesh::Part* ioPart = stk::io::getPart(bulk->mesh_meta_data(), partName);

  EXPECT_TRUE(nullptr != metaPart) << "Could not find meta data part: " << partName;
  EXPECT_TRUE(nullptr !=   ioPart) << "Could not find io part: " << partName;
}

}

