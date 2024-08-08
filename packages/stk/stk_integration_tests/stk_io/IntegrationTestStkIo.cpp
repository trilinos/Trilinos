#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include "stk_util/parallel/Parallel.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "stk_unit_test_utils/ioUtils.hpp"

namespace
{
bool file_exists(const std::string& name)
{
  std::ifstream f(name.c_str());
  return f.good();
}

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

  EXPECT_TRUE(nullptr != metaPart) << "Could not find meta data part: " << partName;
}

TEST(StkIo, checkCanonicalNameFromFile)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) != 1) return;

  stk::io::StkMeshIoBroker stkIo;
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(communicator).create();

  std::string meshSpec = stk::unit_test_util::get_option("-mesh", "");
  std::string partName = stk::unit_test_util::get_option("-part", "UNKNOWN");

  if(file_exists(meshSpec)) {
    stk::io::fill_mesh(meshSpec, *bulk, stkIo);

    stk::mesh::Part* metaPart = bulk->mesh_meta_data().get_part(partName);

    EXPECT_TRUE(nullptr != metaPart) << "Could not find meta data part: " << partName;

    if(nullptr != metaPart) {
      std::cout << "Part name for '" << partName << "' : " << metaPart->name() << std::endl;
    }
  }
}

}

