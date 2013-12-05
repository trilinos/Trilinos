#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>

namespace {

  TEST(StkMeshIoBrokerHowTo, subsetOutputDatabase)
  {
    std::string resultsFilename = "subsetted.results";
    MPI_Comm communicator = MPI_COMM_WORLD;

    size_t num_elems_per_edge = 9;  
    {
      //-BEGIN
      // ============================================================
      // INITIALIZATION
      std::string s_elems_per_edge = Ioss::Utils::to_string(num_elems_per_edge);

      //+ Create a generated mesh containg hexes and shells.
      std::string input_filename = s_elems_per_edge + "x" +
                                   s_elems_per_edge + "x" +
                                   s_elems_per_edge + "|shell:xyzXYZ";

      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.open_mesh_database(input_filename, "generated",
			       stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stkIo.populate_bulk_data();

      stk::mesh::MetaData &meta_data = stkIo.meta_data();
      const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();

      // ============================================================
      //+ EXAMPLE
      //+ Create a selector containing just the shell parts.
      stk::mesh::Selector shell_subset;
      for (size_t i=0; i < all_parts.size(); i++) {
	const stk::mesh::Part *part = all_parts[i];
	stk::topology topo = part->topology();
	if (topo == stk::topology::SHELL_QUAD_4) {
	  shell_subset |= *part;
	}
      }

      // Create the output...
      size_t fh = stkIo.create_output_mesh(resultsFilename,
					   stk::io::WRITE_RESULTS);

      //+ Specify that only the subset of parts selected by the
      //+ "shell_subset" selector will be on the output database.
      stkIo.set_subset_selector(fh, shell_subset);
      stkIo.write_output_mesh(fh);
      // Verification omitted...
      //-END
    }

    {
      // ==================================================
      // VERIFICATION
      // Verify output mesh has correct number of nodes and elements.
      // Note that the output mesh will contain all element blocks;
      // however, the non-shell element block will have zero elements.
      // This is due to the subset_selector subsetting the entities and
      // not the parts...
      Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", resultsFilename,
							 Ioss::READ_MODEL, communicator);
      Ioss::Region ioRegion(iossDb);


      // The output model should consist of the elements and nodes in the 6 shell blocks.
      size_t num_elements = ioRegion.get_property("element_count").get_int();
      size_t num_nodes    = ioRegion.get_property("node_count").get_int();

      // Calculate the expected number of nodes and elements.
      size_t expected_elements = 6 * num_elems_per_edge * num_elems_per_edge ;

      size_t num_nodes_per_edge = num_elems_per_edge+1;
      size_t expected_nodes = 6 * num_nodes_per_edge*num_nodes_per_edge;
      expected_nodes -= 12 * num_nodes_per_edge; // Nodes on each edge were double-counted in previous calculation.
      expected_nodes += 8; // Nodes on each corner were removed in previous calculation; add them back.

      ASSERT_EQ(num_elements, expected_elements);
      ASSERT_EQ(num_nodes, expected_nodes);
    }
    unlink(resultsFilename.c_str());
  }

}
