#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <fieldNameTestUtils.hpp>
#include <restartTestUtils.hpp>

namespace {

TEST(StkMeshIoBrokerHowTo, subsetOutputDatabase)
{
    std::string resultsFilename = "subsetted.results";
    MPI_Comm communicator = MPI_COMM_WORLD;

    {
      std::string input_filename = "9x9x9|shell:xyzXYZ|sideset:xX|nodeset:yY";
      stk::io::StkMeshIoBroker mesh_data(communicator);
      mesh_data.open_mesh_database(input_filename, "generated");
      mesh_data.create_input_mesh();

      // This is done just to define some fields in stk
      mesh_data.populate_bulk_data();

      // Create a selector containing just the shell parts.
      // Due to intimate knowledge of the workings of the generated mesh,
      // these parts will be named "block_2", ... "block_7"
      stk::mesh::MetaData &meta_data = mesh_data.meta_data();
      stk::mesh::Selector shell_subset;
      const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();

      for (size_t i=0; i < all_parts.size(); i++) {
	const stk::mesh::Part *part = all_parts[i];
	stk::topology topo = part->topology();
	if (topo == stk::topology::SHELL_QUAD_4) {
	  shell_subset |= *part;
	}
      }

      // Output...
      size_t fileHandle = mesh_data.create_output_mesh(resultsFilename);
      mesh_data.set_subset_selector(fileHandle, shell_subset);
      mesh_data.write_output_mesh(fileHandle);
    }

    // unlink(resultsFilename.c_str());
}

}
