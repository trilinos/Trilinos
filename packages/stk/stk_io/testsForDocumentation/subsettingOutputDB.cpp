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

    size_t num_elems_per_edge = 9;  // If change this, also change input_filename below.
    {
      std::string input_filename = "9x9x9|shell:xyzXYZ|sideset:xX|nodeset:yY";
      stk::io::StkMeshIoBroker mesh_data(communicator);
      mesh_data.open_mesh_database(input_filename, "generated", stk::io::READ_MESH);
      mesh_data.create_input_mesh();

      mesh_data.populate_bulk_data();

      // Create a selector containing just the shell parts.
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
      size_t fileHandle = mesh_data.create_output_mesh(resultsFilename, stk::io::WRITE_RESULTS);
      mesh_data.set_subset_selector(fileHandle, shell_subset);
      mesh_data.write_output_mesh(fileHandle);
    }

    // ========================================================================
    // Verify output mesh has correct number of nodes and elements.
    // Note that the output mesh will contain all element blocks; however,
    // the non-shell element block will have zero elements.
    // This is due to the subset_selector subsetting the entities and not the parts...
    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", resultsFilename, Ioss::READ_MODEL, communicator);
    Ioss::Region ioRegion(iossDb);

    // The output model should consist of the elements and nodes in the 6 shell blocks.
    size_t num_elements = ioRegion.get_property("element_count").get_int();
    size_t num_nodes    = ioRegion.get_property("node_count").get_int();

    size_t expected_elements = 6 * num_elems_per_edge * num_elems_per_edge ;

    size_t num_nodes_per_edge = num_elems_per_edge+1;
    size_t expected_nodes = 6 * num_nodes_per_edge*num_nodes_per_edge;
    expected_nodes -= 12 * num_nodes_per_edge; // Nodes on each edge were double-counted in previous calculation.
    expected_nodes += 8; // Nodes on each corner were removed in previous calculation; add them back.

    ASSERT_EQ(num_elements, expected_elements);
    ASSERT_EQ(num_nodes, expected_nodes);

    unlink(resultsFilename.c_str());
}

}
