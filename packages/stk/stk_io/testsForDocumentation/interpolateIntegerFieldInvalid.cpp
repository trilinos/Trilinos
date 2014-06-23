#include <gtest/gtest.h>                // for AssertHelper, EXPECT_THROW, etc
#include <unistd.h>                     // for unlink
#include <exception>                    // for exception
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_io/MeshField.hpp"         // for MeshField, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
namespace {

  TEST(StkMeshIoBrokerHowTo, interpolateIntegerFieldInvalid)
  {
    std::string ic_name = "interpolate_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 1) {
      return;
    }

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
