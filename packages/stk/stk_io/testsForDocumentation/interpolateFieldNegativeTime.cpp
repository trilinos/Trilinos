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

  TEST(StkMeshIoBrokerHowTo, interpolateFieldNegativeTime)
  {
    std::string ic_name = "interpolate_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    {
      // ============================================================
      //+ INITIALIZATION:
      //+ Create a mesh with the nodal field "temp" for 3 timesteps.
      //+ The times on the database are -2.0, -1.0, and 0.0
      //+ The value of the field at each node is 0.0 at time 0.0,
      //+ -1.0 at time -1.0, and -2.0 at time -2.0
      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:8x8x8|nodeset:xyz";
      stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature", 1);
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
      for (int i=-2; i <= 0; i++) {
	double time = i;

	for(size_t inode=0; inode<nodes.size(); inode++) {
	  double *fieldData =
	    stk::mesh::field_data(temperature, nodes[inode]);
	  *fieldData = time;
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
      //+ The input mesh database has 3 timesteps with times -2.0, -1.0, 0.0.
      //+ The value of the field "temp" is equal to the time
      //+ Read the "temp" value at times -2.0 to 0.0 with an interval
      //+ of 0.1 (-2.0, -1.9, -1.8, ..., 0.0) and verify that 
      //+ the field contains the correct interpolated value.
      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());

      stkIo.populate_bulk_data();

      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);
      
      // The name of the field on the database is "temp"
      stkIo.add_input_field(stk::io::MeshField(temperature, "temp",
					       stk::io::MeshField::LINEAR_INTERPOLATION));

      for (int i=-20; i <= 0; i++) {
	double time = i/10.0;
	//+ Read the field values from the database and verify that they
	//+ are interpolated correctly.
	stkIo.read_defined_input_fields(time);

	// ============================================================
	//+ VERIFICATION
	// The value of the "temperature" field at all nodes should be 'time'
	for(size_t j=0; j<nodes.size(); j++) {
	  double *fieldData = stk::mesh::field_data(temperature, nodes[j]);
	  EXPECT_DOUBLE_EQ(time, *fieldData);
	}
      }
      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(ic_name.c_str());
  }
}
