// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <iomanip>                      // for operator<<
#include <stk_io/InputFile.hpp>         // for InputFile, etc
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

  TEST(StkMeshIoBrokerHowTo, interpolateFieldCyclic)
  {
    std::string ic_name = "interpolate_field_example.e";
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 1) {
      return;
    }

    {
      // ============================================================
      //+ INITIALIZATION:
      //+ Create a mesh with the nodal field "temp" for 3 timesteps.
      //+ The value of the field at each node is 0.0 at time 0.0,
      //+ 10.0 at time 10.0, and 20.0 at time 20.0
      stk::io::StkMeshIoBroker stkIo(communicator);

      const std::string generatedFileName = "generated:1x1x1|nodeset:xyz";
      stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature =
	stkIo.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field_on_mesh(temperature, stkIo.meta_data().universal_part(),
                                   (stk::mesh::FieldTraits<stk::mesh::Field<double> >::data_type*) nullptr);
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
	double time = i * 10.0;

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
      //+ The input mesh database has 3 timesteps with times 0.0, 10.0, 20.0,
      //+ The value of the field "temp" is equal to the time
      //+ Read the "temp" value at times 0.0 to 10.0 with an interval
      //+ of 0.25 (0.0, 0.25, 0.50, 0.75, ..., 10.0)
      //+ The mapping from analysis time (0.0 to 10.0) to database
      //+ time will be reverse cyclic and scaled.
      //+
      //+ The parameters are:
      //+ * period = 2.0
      //+ * scale = 10.0
      //+ * offset = 0.0
      //+ * cycle type = REVERSING
      //+
      //+ Analysis Time and DB_Time:
      //+ 0  1  2  3   4   5   6   7   8   9  10
      //+ 0 10 20 10   0  10  20  10   0  10  20
      //+ 

      stk::io::StkMeshIoBroker stkIo(communicator);
      size_t idx = stkIo.add_mesh_database(ic_name, stk::io::READ_MESH);
      stkIo.create_input_mesh();

      stk::mesh::Field<double> &temperature = stkIo.
	meta_data().declare_field<stk::mesh::Field<double> >
	(stk::topology::NODE_RANK, "temperature", 1);
      stk::mesh::put_field_on_mesh(temperature, stkIo.meta_data().universal_part(),
                                   (stk::mesh::FieldTraits<stk::mesh::Field<double> >::data_type*) nullptr);

      stkIo.populate_bulk_data();

      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);
      
      // The name of the field on the database is "temp"
      stkIo.add_input_field(stk::io::MeshField(temperature, "temp",
					       stk::io::MeshField::LINEAR_INTERPOLATION));

      //+ Set the periodic parameters on the input mesh...
      double period_length = 2.0;
      double startup = 0.0;
      double scale = 10.0;
      stkIo.get_mesh_database(idx)
	.set_periodic_time(period_length, startup, stk::io::InputFile::REVERSING)
	.set_scale_time(scale)
	.set_start_time(0.0).set_offset_time(0.0).set_stop_time(999.0); // These are optional
      double delta_time = 0.25;
      double time = 0.0;
      double expected = 0.0;
      double exp_inc = 10.0 * delta_time;
      
      while (time <= 10.0) {

	//+ Read the field values from the database and verify that they
	//+ are interpolated correctly.
	stkIo.read_defined_input_fields(time);

	// ============================================================
	//+ VERIFICATION
	// The value of the "temperature" field at all nodes should be 'expected'
	for(size_t i=0; i<nodes.size(); i++) {
	  double *fieldData = stk::mesh::field_data(temperature, nodes[i]);
	  EXPECT_DOUBLE_EQ(expected, *fieldData);
	}
	time += delta_time;
	expected += exp_inc;
	if (expected >= 20.0 || expected <= 0.0) {
	  exp_inc = -exp_inc;
	}
      }
      //-END      
    }
    // ============================================================
    // Cleanup
    unlink(ic_name.c_str());
  }
}
