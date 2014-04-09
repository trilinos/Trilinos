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
#include "Ioss_Field.h"                 // for Field, Field::RoleType, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/IossBridge.hpp"        // for get_field_role
#include "stk_io/MeshField.hpp"         // for MeshField, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/baseImpl/FieldRepository.hpp"  // for FieldVector
#include "stk_topology/topology.hpp"    // for topology, etc
namespace {

  TEST(StkMeshIoBrokerHowTo, restartInterpolatedField)
  {
    std::string rs_name = "restart_interpolate_field.rs";
    std::string ic_name = "interpolate_field.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    // ============================================================
    //+ INITIALIZATION:
    {
      //+ Create a "restart database" with several nodal and element fields,
      //+ and some timesteps...
      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:8x8x8|shell:XYZ|"
	"nodeset:xyz|times:3|variables:nodal,4,element,3,nodeset,2";
      stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.create_input_mesh();
      stkIo.add_all_mesh_fields_as_input_fields();
      stkIo.populate_bulk_data();

      size_t fh = stkIo.create_output_mesh(rs_name, stk::io::WRITE_RESTART);

      // Add all fields as restart fields...
      const stk::mesh::FieldVector &fields = stkIo.meta_data().get_fields();
      for (size_t i=0; i < fields.size(); i++) {
	const Ioss::Field::RoleType* role = stk::io::get_field_role(*fields[i]);
	if ( role && *role == Ioss::Field::TRANSIENT ) {
	  stkIo.add_field(fh, *fields[i]);
	}
      }
      
      // Output the field data.
      for (size_t i=0; i < 3; i++) {
	double time = i;
	stkIo.begin_output_step(fh, time);
	stkIo.write_defined_output_fields(fh);
	stkIo.end_output_step(fh);
      }
    }
    
    {
      //+ Create an "initial condition database" with the nodal field
      //+ "temp" for 10 timesteps - 0.0, 1.0, ..., 9.0.
      //+ The value of the field at each node is the 'time' value.
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
      for (size_t i=0; i < 10; i++) {
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
      //-BEGIN
      // ============================================================
      //+ EXAMPLE: 
      //+ The restart mesh database has 3 timesteps with times 0.0, 1.0, 2.0,
      //+ and several fields.
      //+
      //+ The initial condition database has 10 timesteps with times
      //+ 0.0, 1.0, ..., 9.0 and a nodal variable "temp"
      //+ The value of the field "temp" is equal to the time
      //+
      //+ The example will read the restart database at time 1.0
      //+ and then simulate continuing the analysis at that time
      //+ reading the initial condition data from the other database
      //+ interpolating this data.
      stk::io::StkMeshIoBroker stkIo(communicator);
      size_t ic = stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      size_t rs = stkIo.add_mesh_database(rs_name, stk::io::READ_RESTART);

      //+ "Restart" the calculation...
      double time = 1.0;
      stkIo.set_active_mesh(rs); 
      stkIo.create_input_mesh();

      stkIo.add_all_mesh_fields_as_input_fields();

      stk::mesh::Field<double> &temperature = stkIo.meta_data().
	declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field(temperature, stkIo.meta_data().universal_part());

      // The name of the field on the initial condition database is "temp"
      stkIo.add_input_field(ic, stk::io::MeshField(temperature, "temp",
						   stk::io::MeshField::LINEAR_INTERPOLATION));
      stkIo.populate_bulk_data();

      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

      //+ Read restart data
      stkIo.read_defined_input_fields(time);

      //+ Switch active mesh to "initial condition" database
      stkIo.set_active_mesh(ic); 

      double delta_time = 1.0 / 4.0;
      while (time <= 9.0) {
	//+ Read the field values from the database and verify that they
	//+ are interpolated correctly.
	stkIo.read_defined_input_fields(time);

	// ============================================================
	//+ VERIFICATION
	// The value of the "temperature" field at all nodes should be 'time'
	for(size_t i=0; i<nodes.size(); i++) {
	  double *fieldDataForNode = stk::mesh::field_data(temperature, nodes[i]);
	  EXPECT_DOUBLE_EQ(time, *fieldDataForNode);
	}
	time += delta_time;
      }
      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(rs_name.c_str());
    unlink(ic_name.c_str());
  }
}
