#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetEntities.hpp>
namespace {

  TEST(StkMeshIoBrokerHowTo, interpolateIntegerFieldInvalid)
  {
    std::string ic_name = "interpolate_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    {
      //-BEGIN
      // ============================================================
      //+ EXAMPLE: 
      //+ Interpolated fields cannot be of type integer.
      //+ An exception will be thrown if you try to register an
      //+ integer interpolated field.

      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:8x8x8|nodeset:xyz";
      stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<int> &integer_field = stkIo.meta_data().
	declare_field<stk::mesh::Field<int> >(stk::topology::NODE_RANK, "int_field", 1);
      stk::mesh::put_field(integer_field, stkIo.meta_data().universal_part());
      stkIo.populate_bulk_data();

      EXPECT_THROW(stkIo.add_input_field(stk::io::MeshField(integer_field, "int_field",
			 stk::io::MeshField::LINEAR_INTERPOLATION)), std::exception);
      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(ic_name.c_str());
  }
}
