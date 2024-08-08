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
#include <string>                       // for string, operator+
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, NodeBlockContainer
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, FieldBase
#include "stk_mesh/base/FieldState.hpp"  // for FieldState::StateN, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
namespace Ioss { class DatabaseIO; }

namespace {

TEST(StkMeshIoBrokerHowTo, writeResultsWithMultistateField)
{
  std::string resultsFilename = "output.results";
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  //-BEGIN
  const std::string fieldName = "disp";
  const std::string np1Name = fieldName+"NP1";
  const std::string nName   = fieldName+"N";
  const std::string nm1Name = fieldName+"Nm1";
  {
    // ============================================================
    //+ INITIALIZATION
    const std::string exodusFileName = "generated:1x1x8";
    stk::io::StkMeshIoBroker stkIo(communicator);
    size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    //+ Declare a three-state field
    //+ NOTE: Fields must be declared before "populate_bulk_data()" is called
    //+       since it commits the meta data.
    stk::mesh::Field<double> &field =
        stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, fieldName, 3);
    stk::mesh::put_field_on_mesh(field, stkIo.meta_data().universal_part(), nullptr);

    stkIo.populate_bulk_data();

    size_t fh =
        stkIo.create_output_mesh(resultsFilename, stk::io::WRITE_RESULTS);

    // ============================================================
    //+ EXAMPLE
    //+ Output each state of the multi-state field individually to results file
    stk::mesh::FieldBase *statedFieldNp1 = field.field_state(stk::mesh::StateNP1);
    stk::mesh::FieldBase *statedFieldN   = field.field_state(stk::mesh::StateN);
    stk::mesh::FieldBase *statedFieldNm1 = field.field_state(stk::mesh::StateNM1);

    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

    stkIo.add_field(fh, *statedFieldNp1, np1Name);
    stkIo.add_field(fh, *statedFieldN,   nName);
    stkIo.add_field(fh, *statedFieldNm1, nm1Name);

    // Iterate the application's execute loop five times and output
    // field data each iteration.
    for (int step=0; step < 5; step++) {
      double time = step;

      // Application execution...
      // Generate field data... (details omitted)
      //-END
      double value = 10.0 * time;
      for(size_t i=0; i<nodes.size(); i++) {
        double *np1_data =
            static_cast<double*>(stk::mesh::field_data(*statedFieldNp1, nodes[i]));
        *np1_data = value;
        double *n_data   =
            static_cast<double*>(stk::mesh::field_data(*statedFieldN,   nodes[i]));
        *n_data   = value + 0.1;
        double *nm1_data =
            static_cast<double*>(stk::mesh::field_data(*statedFieldNm1, nodes[i]));
        *nm1_data = value + 0.2;
      }
      //-BEGIN
      //+ Results output...
      stkIo.begin_output_step(fh, time);
      stkIo.write_defined_output_fields(fh);
      stkIo.end_output_step(fh);
    }
    //-END
  }
  //-END

  // ============================================================
  //+ VERIFICATION
  {
    Ioss::DatabaseIO *resultsDb = Ioss::IOFactory::create("exodus", resultsFilename,
                                                          Ioss::READ_MODEL, communicator);
    Ioss::Region results(resultsDb);
    // Should be 5 steps on database...
    EXPECT_EQ(results.get_property("state_count").get_int(), 5);
    // Should be 3 nodal fields on database named "disp.NP1", "disp.N", "disp.NM1";
    Ioss::NodeBlock *nb = results.get_node_blocks()[0];
    EXPECT_EQ(3u, nb->field_count(Ioss::Field::TRANSIENT));
    EXPECT_TRUE(nb->field_exists(np1Name));
    EXPECT_TRUE(nb->field_exists(nName));
    EXPECT_TRUE(nb->field_exists(nm1Name));

    // Iterate each step and verify that the correct data was written.
    for (size_t step=0; step < 5; step++) {
      double time = step;

      double db_time = results.begin_state(step+1);
      EXPECT_EQ(time, db_time);

      std::vector<double> field_data_np1;
      std::vector<double> field_data_n;
      std::vector<double> field_data_nm1;
      nb->get_field_data(np1Name, field_data_np1);
      nb->get_field_data(nName,   field_data_n);
      nb->get_field_data(nm1Name, field_data_nm1);

      EXPECT_EQ(field_data_np1.size(), field_data_n.size());
      EXPECT_EQ(field_data_np1.size(), field_data_nm1.size());

      double expected = 10.0 * time;
      for (size_t node = 0; node < field_data_np1.size(); node++) {
        EXPECT_EQ(field_data_np1[node], expected);
        EXPECT_EQ(field_data_n[node],   expected+0.1);
        EXPECT_EQ(field_data_nm1[node], expected+0.2);
      }
      results.end_state(step+1);
    }
  }

  unlink(resultsFilename.c_str());
}
}
