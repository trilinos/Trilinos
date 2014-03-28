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

  TEST(StkMeshIoBrokerHowTo, handleMissingFieldOnRead)
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

      const std::string generatedFileName = "generated:8x8x8|nodeset:xyz";
      size_t index = stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.set_active_mesh(index); // Optional if only a single input database
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,"temperature",1);
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
	  double *fieldDataForNode =
	    stk::mesh::field_data(temperature, nodes[inode]);
	  *fieldDataForNode = time;
	}

	stkIo.begin_output_step(fh, time);
	stkIo.write_defined_output_fields(fh);
	stkIo.end_output_step(fh);
      }
    }

    {
      // ============================================================
      //+ EXAMPLE: 
      //+ Demonstrate what happens when application requests the
      //+ reading of a field that does not exist on the input
      //+ mesh database.  The nodal field "displacement" is
      //+ requested for input from the database field "disp" which
      //+ does not exist.
      stk::io::StkMeshIoBroker stkIo(communicator);
      size_t index = stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      stkIo.set_active_mesh(index);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,"temperature",1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());

      stk::mesh::Field<double> &displacement = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,"displacement",1);
      stk::mesh::put_field(displacement, stkIo.meta_data().universal_part());
      stkIo.populate_bulk_data();

      // The name of the field on the database is "temp"
      // This field does exist and should be read correctly
      stkIo.add_input_field(stk::io::MeshField(temperature, "temp"));

      //+ The name of the field on the database is "disp"
      //+ This field does not exist and will not be found.
      stkIo.add_input_field(stk::io::MeshField(displacement, "disp"));

      //-BEGIN
      //+ If read the fields, but don't pass in the 'missing_fields'
      //+ vector, the code will print an error message and throw an
      //+ exception if it doesn't find all of the requested fields.
      EXPECT_THROW(stkIo.read_defined_input_fields(2.0), std::exception);

      //+ If code throws due to missing field(s), it will NOT read
      //+ even the fields that exist.
      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(ic_name.c_str());
  }
}
