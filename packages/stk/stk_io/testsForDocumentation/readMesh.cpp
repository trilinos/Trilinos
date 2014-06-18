#include <gtest/gtest.h>                // for AssertHelper, EXPECT_NE, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stddef.h>                     // for size_t, NULL
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_util/parallel/Parallel.hpp"
namespace stk { namespace mesh { class Part; } }
namespace {

  TEST(StkMeshIoBrokerHowTo, readMesh)
  {
    std::string mesh_name = "input_mesh_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 1) {
      return;
    }

    {
      // ============================================================
      //+ INITIALIZATION:
      //+ Create a basic mesh with a hex block, 3 shell blocks, 3 nodesets, and 3 sidesets.
      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:8x8x8|shell:xyz|nodeset:xyz|sideset:XYZ";
      size_t index = stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.set_active_mesh(index);

      stkIo.create_input_mesh();
      stkIo.populate_bulk_data();

      size_t fh = stkIo.create_output_mesh(mesh_name, stk::io::WRITE_RESULTS);
      stkIo.write_output_mesh(fh);
    }

    {
      //-BEGIN
      // ============================================================
      //+ EXAMPLE: 
      //+ Read mesh data from the specified file.
      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.add_mesh_database(mesh_name, stk::io::READ_MESH);

      //+ Creates meta data; creates parts 
      stkIo.create_input_mesh();

      //+ Any modifications to the meta data must be done here.
      //+ This includes declaring fields.
      
      //+ Commit the meta data and create the bulk data.  
      //+ populate the bulk data with data from the mesh file.
      stkIo.populate_bulk_data();

      // ============================================================
      //+ VERIFICATION
      //+ There should be:
      //+ 4 parts corresponding to the 1 hex block and 3 shell blocks
      stk::mesh::MetaData &meta = stkIo.meta_data();
      stk::mesh::Part *invalid = NULL;
      EXPECT_NE(invalid, meta.get_part("block_1"));
      EXPECT_NE(invalid, meta.get_part("block_2"));
      EXPECT_NE(invalid, meta.get_part("block_3"));
      EXPECT_NE(invalid, meta.get_part("block_4"));

      //+ 3 parts corresponding to the 3 nodesets.
      EXPECT_NE(invalid, meta.get_part("nodelist_1"));
      EXPECT_NE(invalid, meta.get_part("nodelist_2"));
      EXPECT_NE(invalid, meta.get_part("nodelist_3"));

      //+ 3 parts corresponding to the 3 sidesets.
      EXPECT_NE(invalid, meta.get_part("surface_1"));
      EXPECT_NE(invalid, meta.get_part("surface_2"));
      EXPECT_NE(invalid, meta.get_part("surface_3"));

      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(mesh_name.c_str());
  }
}
