#include <gtest/gtest.h>
#include <stk_unit_test_utils/getOption.h>
#include <unistd.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_unit_test_utils/CommandLineArgs.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>

namespace
{

TEST(StkIoHowTo, WriteMesh)
{
  std::string filename = "output.exo";
  {
    std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
    stk::io::fill_mesh("generated:1x1x4", *bulk);

    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(bulk);
    size_t outputFileIndex = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(outputFileIndex);
    stkIo.write_defined_output_fields(outputFileIndex);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
    stk::io::fill_mesh(filename, *bulk);

    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(*bulk, entityCounts);
    EXPECT_EQ(4u, entityCounts[stk::topology::ELEM_RANK]);
  }

  unlink(filename.c_str());
}

TEST(StkIoHowTo, generateMeshWith64BitIds)
{
  std::string meshSpec = stk::unit_test_util::get_option("-i", "1x1x4");
  std::string fullMeshSpec = "generated:"+meshSpec;

  std::string filename = "output.exo";
  stk::io::StkMeshIoBroker inputBroker;

  //+ Set properties to ensure that 64-bit integers will be used
  inputBroker.property_add(Ioss::Property("INTEGER_SIZE_API" , 8));
  inputBroker.property_add(Ioss::Property("INTEGER_SIZE_DB" , 8));
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::io::fill_mesh_preexisting(inputBroker, fullMeshSpec, *bulk);

  stk::io::write_mesh_with_large_ids_and_fields(filename, *bulk);

  if (stk::unit_test_util::GlobalCommandLineArguments::self().get_argc() == 0) {
    unlink(filename.c_str());
  }
}

}
