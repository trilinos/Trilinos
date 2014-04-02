#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Types.hpp>
#include <Ioss_SubSystem.h>

namespace {

  TEST(StkMeshIoBrokerHowTo, writeResults)
  {
    std::string mesh_name = "input_mesh_example.e";
    std::string results_name = "output.results";
    MPI_Comm communicator = MPI_COMM_WORLD;

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
      // ============================================================
      //+ EXAMPLE: 
      //+ Read mesh data from the specified file.
      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.add_mesh_database(mesh_name, stk::io::READ_MESH);

      //+ Creates meta data; creates parts 
      stkIo.create_input_mesh();

      //+ Declare a field
      //+ NOTE: Fields must be declared before "populate_bulk_data()" is called
      //+       since it commits the meta data.
      const std::string fieldName = "disp";
      stk::mesh::Field<double> &field =
	stkIo.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,
								   fieldName, 1);
      stk::mesh::put_field(field, stkIo.meta_data().universal_part());

      //+ commit the meta data and create the bulk data.  
      //+ populate the bulk data with data from the mesh file.
      stkIo.populate_bulk_data();

      // ============================================================
      //+ Create results file. By default, all parts created from the input
      //+ mesh will be written to the results output file.
      size_t fh = stkIo.create_output_mesh(results_name, stk::io::WRITE_RESULTS);

      //-BEGIN
      //+ The field 'fieldName' will be output to the results file with the name 'alternateFieldName'
      std::string alternateFieldName("displacement");
      stkIo.add_field(fh, field, alternateFieldName);
      //-END      

      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

      // Iterate the application's execute loop five times and output
      // field data each iteration.
      for (int step=0; step < 5; step++) {
	double time = step;

	// Application execution...
	double value = 10.0 * time;
	for(size_t i=0; i<nodes.size(); i++) {
	  double *node_data = stk::mesh::field_data(field, nodes[i]);
	  *node_data = value;
	}

	//+ Output the field data calculated by the application.
	stkIo.begin_output_step(fh, time);
	stkIo.write_defined_output_fields(fh);
        stkIo.end_output_step(fh);
      }

      // ============================================================
      //+ VERIFICATION
      Ioss::Region *io_region = stkIo.get_output_io_region(fh).get();
      Ioss::NodeBlock *nb = io_region->get_node_blocks()[0];
      ASSERT_TRUE(nb->field_exists(alternateFieldName));
      ASSERT_FALSE(nb->field_exists(fieldName));
    }

    // ============================================================
    // Cleanup
    unlink(results_name.c_str());
  }
}

