#include <gtest/gtest.h>                // for AssertHelper, etc
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
#include "stk_util/parallel/Parallel.hpp"
namespace {

  TEST(StkMeshIoBrokerHowTo, interpolateOutsideRange)
  {
    std::string ic_name = "interpolate_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;
    int num_procs = stk::parallel_machine_size(communicator);
    if (num_procs != 1) {
      return;
    }

    {
      // ============================================================
      //+ INITIALIZATION:
      //+ Create a mesh with the nodal field "temp" for 2 timesteps
      //+ with times 1.0 and 2.0.
      //+ The value of the field at each node is equal to the 'time'
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
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);
    
      // Add two steps to the database
      for (size_t i=1; i <= 2; i++) {
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
      //+ The input mesh database has 2 timesteps with time 1.0 and 2.0
      //+ The value of the field "temp" is equal to the time
      //+ Read the "temp" value at times 0.0 to 3.0 with an interval
      //+ of 0.1 (0.0, 0.1, 0.2, 0.3, ..., 2.0).
      //+
      //+ The times 0.0 to 1.0 and 2.0 to 3.0 are outside
      //+ the range of the mesh database so no interpolation
      //+ or extrapolation will occur -- the field values
      //+ will be set to the values at the nearest time.
      //+
      //+ Verify that the values from times 0.0 to 1.0
      //+ are equal to 1.0 and that the values from 2.0 to 3.0
      //+ are equal to 2.0.
      //+ The field values from 1.0 to 2.0 will be interpolated
      //+
      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,"temperature",1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());

      stkIo.populate_bulk_data();

      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

      // The name of the field on the database is "temp"
      stkIo.add_input_field(stk::io::MeshField(temperature, "temp",
					       stk::io::MeshField::LINEAR_INTERPOLATION));

      for (size_t i=0; i < 21; i++) {
	double time = i/10.0;
	//+ Read the field values from the database and verify that they
	//+ are interpolated correctly.
	stkIo.read_defined_input_fields(time);

	// ============================================================
	//+ VERIFICATION

	double expected_value = time;
	if (time <= 1.0)
	  expected_value = 1.0;
	if (time >= 2.0)
	  expected_value = 2.0;
		
	for(size_t j=0; j<nodes.size(); j++) {
	  double *fieldData = stk::mesh::field_data(temperature, nodes[j]);
	  EXPECT_DOUBLE_EQ(expected_value, *fieldData);
	}
      }
      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(ic_name.c_str());
  }
}
