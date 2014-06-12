#include <gtest/gtest.h>                // for AssertHelper, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <iomanip>                      // for operator<<
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/MeshField.hpp"         // for MeshField, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_topology/topology.hpp"    // for topology, etc
namespace {

  TEST(StkMeshIoBrokerHowTo, interpolateNodalField)
  {
    std::string ic_name = "interpolate_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    {
      // ============================================================
      //+ INITIALIZATION:
      //+ Create a mesh with the nodal field "temp" for 3 timesteps.
      //+ The value of the field at each node is 0.0 at time 0.0,
      //+ 1.0 at time 1.0, and 2.0 at time 2.0
      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:8x8x8|nodeset:xyz";
      stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
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
      //+ The input mesh database has 3 timesteps with times 0.0, 1.0, 2.0,
      //+ The value of the field "temp" is equal to the time
      //+ Read the "temp" value at times 0.0 to 2.0 with an interval
      //+ of 0.1 (0.0, 0.1, 0.2, 0.3, ..., 2.0) and verify that 
      //+ the field contains the correct interpolated value.
      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,"temperature",1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());

      stkIo.populate_bulk_data();

      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);
      
      //+ Specify that the field data are to be linear interpolated.
      stkIo.add_input_field(stk::io::MeshField(temperature, "temp",
					       stk::io::MeshField::LINEAR_INTERPOLATION));/*@\label{io:input:interpolate}*/

      //+ If the same stk field (temperature) is added more than once,
      //+ the first database name and settings will be used.  For example,
      //+ the add_input_field below will be ignored with no error or warning.
      stkIo.add_input_field(stk::io::MeshField(temperature, "temp-again",
					       stk::io::MeshField::LINEAR_INTERPOLATION));

      for (size_t i=0; i < 21; i++) {
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
