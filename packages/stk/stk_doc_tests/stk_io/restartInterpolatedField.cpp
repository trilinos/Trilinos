// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
#include "stk_mesh/base/Types.hpp"      // for FieldVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
namespace {

TEST(StkMeshIoBrokerHowTo, restartInterpolatedField)
{
  std::string rs_name = "restart_interpolate_field.rs";
  std::string ic_name = "interpolate_field.e";
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

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

    stk::mesh::Field<double> &temperature =
        stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, "temperature", 1);
    stk::mesh::put_field_on_mesh(temperature, stkIo.meta_data().universal_part(), nullptr);
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

    stk::mesh::Field<double> &temperature =
        stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, "temperature", 1);
    stk::mesh::put_field_on_mesh(temperature, stkIo.meta_data().universal_part(), nullptr);

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
