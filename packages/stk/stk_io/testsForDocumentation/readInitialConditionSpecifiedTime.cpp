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

  TEST(StkMeshIoBrokerHowTo, readInitialConditionSpecifiedTime)
  {
    std::string ic_name = "input_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    {
      // ============================================================
      //+ INITIALIZATION:
      //+ Create a mesh with the 2 nodal fields "temp" and "flux"
      //+ for 3 timesteps.
      //+ The value of the fields at each node is 0.0 at time 0.0,
      //+ 1.0 at time 1.0, and 2.0 at time 2.0
      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:8x8x8|nodeset:xyz";
      stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());

      stk::mesh::Field<double> &heat_flux = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "heat_flux", 1);
      stk::mesh::put_field(heat_flux, stkIo.meta_data().universal_part());
      stkIo.populate_bulk_data();

      size_t fh = stkIo.create_output_mesh(ic_name, stk::io::WRITE_RESULTS);

      stkIo.add_field(fh, temperature, "temp");
      stkIo.add_field(fh, heat_flux, "flux");
	
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

	  fieldDataForNode = stk::mesh::field_data(heat_flux, nodes[inode]);
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
      //+ Register the reading of database fields "temp" and "flux" to
      //+ populate the stk nodal fields "temperature" and "heat_flux"
      //+ for use as initial conditionss.
      //+ Specify that the "temp" field should be read from database
      //+ time 2.0 no matter what time is specified in the read_defined_input_fields
      //+ call.
      //+ The "flux" field will be read at the database time corresponding
      //+ to the analysis time passed in to read_defined_input_fields.

      stk::io::StkMeshIoBroker stkIo(communicator);
      size_t index = stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      stkIo.set_active_mesh(index);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());

      stk::mesh::Field<double> &heat_flux = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "heat_flux", 1);
      stk::mesh::put_field(heat_flux, stkIo.meta_data().universal_part());
      stkIo.populate_bulk_data();

      // The name of the field on the database is "temp"
      stk::io::MeshField temp_field(temperature, "temp", stk::io::MeshField::CLOSEST);
      temp_field.set_read_time(2.0);/*@\label{io:input:specifiedtime}*/
      stkIo.add_input_field(temp_field);

      // The name of the field on the database is "flux"
      stk::io::MeshField flux_field(heat_flux, "flux", stk::io::MeshField::CLOSEST);
      stkIo.add_input_field(flux_field);

      //+ Read the field values from the database at time 1.0
      //+ The value of "flux" will be the values from database time 1.0
      //+ However, the value of "temp" will be the values from database time 2.0
      stkIo.read_defined_input_fields(1.0);

      // ============================================================
      //+ VERIFICATION
      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK,
			      nodes);

      //+ The value of the "temperature" field at all nodes should be 2.0
      for(size_t i=0; i<nodes.size(); i++) {
	double *fieldDataForNode = stk::mesh::field_data(temperature, nodes[i]);
	EXPECT_DOUBLE_EQ(2.0, *fieldDataForNode);
      }

      //+ The value of the "heat_flux" field at all nodes should be 1.0
      for(size_t i=0; i<nodes.size(); i++) {
	double *fieldDataForNode = stk::mesh::field_data(heat_flux, nodes[i]);
	EXPECT_DOUBLE_EQ(1.0, *fieldDataForNode);
      }
      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(ic_name.c_str());
  }
}
