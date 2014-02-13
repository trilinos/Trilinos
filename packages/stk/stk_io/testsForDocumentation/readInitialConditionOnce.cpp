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

  TEST(StkMeshIoBrokerHowTo, readInitialConditionOnce)
  {
    std::string ic_name = "input_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    {
      // ============================================================
      //+ INITIALIZATION:
      //+ Create a mesh with the nodal field "temp" for 3 timesteps.
      //+ The value of the field at each node is 0.0 at time 0.0,
      //+ 1.0 at time 1.0, and 2.0 at time 2.0
      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:8x8x8";
      size_t index = stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.set_active_mesh(index);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature =
	stkIo.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());
      stkIo.populate_bulk_data();

      size_t fh = stkIo.create_output_mesh(ic_name, stk::io::WRITE_RESULTS);

      //+ The name of the field on the database will be "temp"
      stkIo.add_field(fh, temperature, "temp");
	
      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(),
			      stk::topology::NODE_RANK, nodes);
    
      // Add three steps to the database
      // For each step, the value of the field is the value 'time'
      for (size_t i=0; i < 3; i++) {
	double time = i;

	for(size_t inode=0; inode<nodes.size(); inode++) {
	  double *fieldDataForNode = stk::mesh::field_data(temperature, nodes[inode]);
	  *fieldDataForNode = time;
	}

	stkIo.begin_output_step(fh, time);
	stkIo.write_defined_output_fields(fh);
	stkIo.end_output_step(fh);
      }
    }

    {
      //-BEGIN
      // ============================================================
      //+ EXAMPLE: 
      //+ Read the value of the "temp" field at step 2 and populate
      //+ the nodal field "temperature" for use as an initial condition
      //+ The input field should only be active for one 'read_defined_input_fields'
      //+ call, so verify this by calling the function again at step 3 and
      //+ then verify that the field values are still those read from step 2.
      stk::io::StkMeshIoBroker stkIo(communicator);
      size_t index = stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      stkIo.set_active_mesh(index);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature =
	stkIo.meta_data().declare_field<stk::mesh::Field<double> >
	                                            (stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());
      stkIo.populate_bulk_data();

      //+ The name of the field on the database is "temp"
      stk::io::MeshField input_field(temperature, "temp", stk::io::MeshField::CLOSEST);
      input_field.set_read_once(true);
      stkIo.add_input_field(input_field);

      //+ Read the field values from the database at time 2.0
      //+ Pass in a time of 2.2 to verify that the value returned is
      //+ from the closest step and not interpolated.
      stkIo.read_defined_input_fields(2.2);

      // ============================================================
      //+ VERIFICATION
      //+ The value of the field at all nodes should be 2.0
      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK,
			      nodes);
      for(size_t i=0; i<nodes.size(); i++) {
	double *fieldDataForNode = stk::mesh::field_data(temperature, nodes[i]);
	EXPECT_DOUBLE_EQ(2.0, *fieldDataForNode);
      }

      //+ Call read_defined_input_fields again and verify that the
      //+ input field registration is no longer active after the
      //+ since it was specified to be "only_read_once()"
      stkIo.read_defined_input_fields(3.0);

      // ============================================================
      //+ VERIFICATION
      //+ The value of the field at all nodes should still be 2.0
      for(size_t i=0; i<nodes.size(); i++) {
	double *fieldDataForNode = stk::mesh::field_data(temperature, nodes[i]);
	EXPECT_DOUBLE_EQ(2.0, *fieldDataForNode);
      }

      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(ic_name.c_str());
  }
}
