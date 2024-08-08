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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <string>                       // for string
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, NodeBlockContainer
#include "Ioss_ElementBlock.h"
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, FieldBase
#include "stk_mesh/base/FieldState.hpp"  // for FieldState::StateN, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace Ioss { class DatabaseIO; }
namespace {

TEST(StkMeshIoBrokerHowTo, writeResultsAndRestart)
{
  std::string mesh_name = "input_mesh_example.e";
  std::string results_name = "output.results";
  std::string restart_name = "output.restart";
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
    //-BEGIN
    // ============================================================
    //+ EXAMPLE:
    //+ Read mesh data from the specified file.
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(mesh_name, stk::io::READ_MESH);

    //+ Creates meta data; creates parts
    stkIo.create_input_mesh();

    //+ Declare a three-state field
    //+ NOTE: Fields must be declared before "populate_bulk_data()" is called
    //+       since it commits the meta data.
    const std::string fieldName = "disp";
    stk::mesh::Field<double> &field = stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, fieldName, 3);
    stk::mesh::put_field_on_mesh(field, stkIo.meta_data().universal_part(), nullptr);

    const stk::mesh::Part& block_1 = *stkIo.meta_data().get_part("block_1");
    //+ create a two-state field
    stk::mesh::Field<double> &fooSubset = stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, "fooSubset", 2);
    stk::mesh::put_field_on_mesh(fooSubset, block_1, nullptr);

    //+ commit the meta data and create the bulk data.
    //+ populate the bulk data with data from the mesh file.
    stkIo.populate_bulk_data();

    // ============================================================
    //+ Create results file. By default, all parts created from the input
    //+ mesh will be written to the results output file.
    size_t results_fh = stkIo.create_output_mesh(results_name, stk::io::WRITE_RESULTS);

    //+ Create restart file. By default, all parts created from the input
    //+ mesh will be written to the results output file.
    size_t restart_fh = stkIo.create_output_mesh(restart_name, stk::io::WRITE_RESTART);

    //+ The field will be output to the results file with the default field name.
    //+ Only the newest state will be output.
    stkIo.add_field(results_fh, field);

    //+ Output the field to the restart database also.
    //+ The two newest states will be output.
    stkIo.add_field(restart_fh, field);
    stkIo.add_field(restart_fh, fooSubset);

    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

    stk::mesh::FieldBase *statedFieldNp1 = field.field_state(stk::mesh::StateNP1);
    stk::mesh::FieldBase *statedFieldN   = field.field_state(stk::mesh::StateN);
    stk::mesh::FieldBase *statedFieldNm1 = field.field_state(stk::mesh::StateNM1);

    // Iterate the application's execute loop five times and output
    // field data each iteration.
    for (int step=0; step < 5; step++) {
      double time = step;

      // Application execution...
      double value = 10.0 * time;
      for(size_t i=0; i<nodes.size(); i++) {
        double *np1_data = static_cast<double*>(stk::mesh::field_data(*statedFieldNp1, nodes[i]));
        *np1_data = value;
        double *n_data   = static_cast<double*>(stk::mesh::field_data(*statedFieldN,   nodes[i]));
        *n_data   = value + 0.1;
        double *nm1_data = static_cast<double*>(stk::mesh::field_data(*statedFieldNm1, nodes[i]));
        *nm1_data = value + 0.2;
      }

      //+ Results output...
      stkIo.begin_output_step(results_fh, time);
      stkIo.write_defined_output_fields(results_fh);
      stkIo.end_output_step(results_fh);

      //+ Restart output...
      stkIo.begin_output_step(restart_fh, time);
      stkIo.write_defined_output_fields(restart_fh);
      stkIo.end_output_step(restart_fh);
    }
    //-END
  }

  {
    //Demonstrate reading a restart database after adding an extra multistate field to
    //the io-broker "foo" which we know is missing from the database. We confirm that
    //passing the 'missingFields' argument to stkIo.read_defined_input_fields allows
    //the code to continue without throwing an exception due to not finding the field.
    stk::io::StkMeshIoBroker stkIo(communicator);
    size_t rs = stkIo.add_mesh_database(restart_name, stk::io::READ_RESTART);

    //+ "Restart" the calculation...
    double time = 1.0;
    stkIo.set_active_mesh(rs);
    stkIo.create_input_mesh();

    //create a 3-state field
    stk::mesh::Field<double> &foo =
        stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, "foo", 3);
    stk::mesh::put_field_on_mesh(foo, stkIo.meta_data().universal_part(), nullptr);

    const stk::mesh::Part& block_1 = *stkIo.meta_data().get_part("block_1");
    //create a 3-state field
    stk::mesh::Field<double> &fooSubset = stkIo.meta_data().
        declare_field<double>(stk::topology::NODE_RANK, "fooSubset", 3);
    stk::mesh::put_field_on_mesh(fooSubset, block_1, nullptr);

    //add the new fields to the stk-io-broker, even though we know it isn't present
    //in the restart database. This is to test the 'missing-fields' argument below.
    stk::io::MeshField meshfieldFoo(&foo, "FOO");
    meshfieldFoo.set_single_state(false);
    meshfieldFoo.add_part(stk::topology::NODE_RANK, stkIo.meta_data().universal_part(),
                          stkIo.get_input_ioss_region().get()->get_node_blocks()[0]);
    stkIo.add_input_field(meshfieldFoo);

    stk::io::MeshField meshfieldFooSubset(&fooSubset, "FOOSUBSET");
    meshfieldFooSubset.set_single_state(false);
    const Ioss::ElementBlock& elem_block = *stkIo.get_input_ioss_region().get()->get_element_blocks()[0];
    Ioss::ElementBlock& nonconst_elem_block = const_cast<Ioss::ElementBlock&>(elem_block);
    meshfieldFooSubset.add_part(stk::topology::NODE_RANK, block_1, &nonconst_elem_block);
    meshfieldFooSubset.add_subset(block_1);
    stkIo.add_input_field(meshfieldFooSubset);

    stk::io::set_field_role(foo, Ioss::Field::TRANSIENT);
    stk::io::set_field_role(fooSubset, Ioss::Field::TRANSIENT);

    //stkIo.add_all_mesh_fields_as_input_fields();

    stkIo.populate_bulk_data();

    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

    std::vector<stk::io::MeshField> missingFields;
    //+ Read restart data
    stkIo.read_defined_input_fields(time, &missingFields);

    EXPECT_EQ(3u, missingFields.size());
    const stk::io::MeshField& missingField0 = missingFields[0];
    std::string name = missingField0.db_name();
    std::string statedFieldName = missingField0.field()->name();
    EXPECT_EQ("FOO", name);
    EXPECT_EQ("foo", statedFieldName);
    const stk::io::MeshField& missingField1 = missingFields[1];
    name = missingField1.db_name();
    statedFieldName = missingField1.field()->name();
    EXPECT_EQ("FOO", name);
    EXPECT_EQ("foo_STKFS_N", statedFieldName);
    const stk::io::MeshField& missingField2 = missingFields[2];
    statedFieldName = missingField2.field()->name();
    name = missingField2.db_name();
    EXPECT_EQ("FOOSUBSET", name);
    EXPECT_EQ("fooSubset_STKFS_N", statedFieldName);
  }

  // ============================================================
  //+ VERIFICATION
  {
    Ioss::DatabaseIO *resultsDb = Ioss::IOFactory::create("exodus", results_name,
                                                          Ioss::READ_MODEL, communicator);
    Ioss::Region results(resultsDb);
    // Should be 5 steps on database...
    EXPECT_EQ(results.get_property("state_count").get_int(), 5);
    // Should be 1 nodal field on database named "disp";
    Ioss::NodeBlock *nb = results.get_node_blocks()[0];
    EXPECT_EQ(1u, nb->field_count(Ioss::Field::TRANSIENT));
    EXPECT_TRUE(nb->field_exists("disp"));

    // Iterate each step and verify that the correct data was written.
    for (size_t step=0; step < 5; step++) {
      double time = step;

      double db_time = results.begin_state(step+1);
      EXPECT_EQ(time, db_time);

      std::vector<double> field_data;
      nb->get_field_data("disp", field_data);
      double expected = 10.0 * time;
      for (size_t node = 0; node < field_data.size(); node++) {
        EXPECT_EQ(field_data[node], expected);
      }
      results.end_state(step+1);
    }
  }

  {
    Ioss::DatabaseIO *restartDb = Ioss::IOFactory::create("exodus", restart_name,
                                                          Ioss::READ_MODEL, communicator);
    Ioss::Region restart(restartDb);
    // Should be 5 steps on database...
    EXPECT_EQ(restart.get_property("state_count").get_int(), 5);
    // Should be 2 nodal field on database named "disp" and "disp.N";
    Ioss::NodeBlock *nb = restart.get_node_blocks()[0];
    EXPECT_EQ(3u, nb->field_count(Ioss::Field::TRANSIENT));
    EXPECT_TRUE(nb->field_exists("disp"));
    EXPECT_TRUE(nb->field_exists("disp.N"));

    // Iterate each step and verify that the correct data was written.
    for (size_t step=0; step < 5; step++) {
      double time = step;

      double db_time = restart.begin_state(step+1);
      EXPECT_EQ(time, db_time);

      std::vector<double> field_data_n;
      std::vector<double> field_data_np1;
      nb->get_field_data("disp", field_data_np1);
      nb->get_field_data("disp.N", field_data_n);
      double expected = 10.0 * time;
      for (size_t node = 0; node < field_data_n.size(); node++) {
        EXPECT_EQ(field_data_np1[node], expected);
        EXPECT_EQ(field_data_n[node],   expected+0.1);
      }
      restart.end_state(step+1);
    }
  }

  // ============================================================
  // Cleanup
  unlink(mesh_name.c_str());
  unlink(results_name.c_str());
  unlink(restart_name.c_str());
}
}
