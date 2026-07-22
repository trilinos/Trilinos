
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

#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_ElementBlock.h"          // for ElementBlock
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, ElementBlockContainer
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/FillMesh.hpp"
#include "stk_io/WriteMesh.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/FieldBLAS.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <ostream>                      // for basic_ostream::operator<<
#include <stddef.h>                     // for size_t
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include "stk_util/parallel/ParallelReduce.hpp"
#include <string>                       // for string
#include <unistd.h>                     // for unlink
#include <vector>                       // for vector
namespace Ioss { class DatabaseIO; }

namespace {

TEST(StkMeshIoBrokerHowTo, writeResultsWithDTIO)
{
  std::string mesh_name = "input_mesh_example.e";
  std::string results_name = "output.e";
  MPI_Comm communicator = MPI_COMM_WORLD;

  {
    // ============================================================
    //+ INITIALIZATION:
    //+ Create a basic mesh with a hex block, 3 shell blocks, 3 nodesets, and 3 sidesets.
    std::unique_ptr<stk::mesh::BulkData> mesh = stk::mesh::MeshBuilder(communicator).create();

    const std::string generatedFileName = "generated:2x2x2|shell:xyz|nodeset:xyz|sideset:XYZ";
    stk::io::fill_mesh(generatedFileName, *mesh);
    stk::io::write_mesh(mesh_name, *mesh);
  }

  {
    //-BEGIN
    // ============================================================
    //+ EXAMPLE:
    //+ Read mesh data from the specified file.
    stk::io::StkMeshIoBroker stkIo(communicator);

    //+ Enable automatic IOSS control of dynamic topology IO (DTIO) to create a single group file
    stkIo.enable_dynamic_topology(stk::io::FileOption::USE_DYNAMIC_TOPOLOGY_GROUP_FILE);

    stkIo.add_mesh_database(mesh_name, stk::io::READ_MESH);

    //+ Creates meta data; creates parts
    stkIo.create_input_mesh();

    //+ Declare a field
    //+ NOTE: Fields must be declared before "populate_bulk_data()" is called
    //+       since it commits the meta data.
    const std::string fieldName = "density";
    stk::mesh::Field<double> &field = stkIo.meta_data().declare_field<double>(stk::topology::ELEM_RANK, fieldName, 1);
    stk::mesh::put_field_on_mesh(field, stkIo.meta_data().universal_part(), nullptr);

    //+ commit the meta data and create the bulk data.
    //+ populate the bulk data with data from the mesh file.
    stkIo.populate_bulk_data();

    // ============================================================
    //+ Create results file. By default, all parts created from the input
    //+ mesh will be written to the results output file.
    size_t fh = stkIo.create_output_mesh(results_name, stk::io::WRITE_RESULTS);

    //+ The field will be output to the results file with the default field name.
    stkIo.add_field(fh, field); /*@\label{io:results:add_field}*/

    //+ Iterate the application's execute loop nine times and output
    //+ field data each iteration (mimicing time steps). The last eight
    //+ "time steps" involve a dynamic topology change where an element is
    //+ destroyed and would normally necessitate the creation of a new
    //+ file. DTIO allows the automatic creation of a single database
    //+ that contains all the files
    for (int step=0; step < 9; step++) {
      double time = step;

      //simulate time-varying result field...
      double value = 10.0 * time;
      stk::mesh::field_fill(value, field);

      if(step != 0) {
        // Destroy an element
        stk::mesh::BulkData& bulk = stkIo.bulk_data();

        bulk.modification_begin();
        stk::mesh::EntityId entityId = step;
        stk::mesh::Entity entity = bulk.get_entity(stk::topology::ELEM_RANK, entityId);
        EXPECT_TRUE(bulk.is_valid(entity));
        EXPECT_TRUE(bulk.destroy_entity(entity));
        bulk.modification_end();

        // Notify STKIO that a mesh change occured
        stkIo.set_topology_modification(fh, Ioss::TOPOLOGY_CREATEELEM);
      }
      //+ Output the field data calculated by the application.
      stkIo.begin_output_step(fh, time);
      stkIo.write_defined_output_fields(fh);
      stkIo.end_output_step(fh);
    }
    //-END
  }
  // ============================================================
  //+ VERIFICATION
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator)
                                                 .set_spatial_dimension(3).create();
    stk::io::StkMeshIoBroker stkIo;

    // Make sure IOSS retains empty element blocks ... final output has no HEX elements
    stkIo.property_add(Ioss::Property("RETAIN_EMPTY_BLOCKS", "YES"));

    stk::io::fill_mesh_preexisting(stkIo, results_name, *bulkPtr);

    // Should be 9 mesh groups in database...one per output time step
    constexpr int goldNumMeshGroups = 9;
    EXPECT_EQ(goldNumMeshGroups, stkIo.num_mesh_groups());

    for(int i=0; i<goldNumMeshGroups; i++) {
      // Load each individual mesh group
      // First one is actually a no-op since group 0 is loaded by default
      EXPECT_TRUE(stkIo.load_mesh_group(i));

      // One time step per mesh group was written
      constexpr int goldNumTimeSteps = 1;
      EXPECT_EQ(goldNumTimeSteps, stkIo.get_num_time_steps());

      const auto& region = stkIo.get_input_ioss_region();
      const auto& eBlocks = region->get_element_blocks();

      // 1 HEX and 3 SHELL blocks
      EXPECT_EQ(4u, eBlocks.size());

      // Look at HEX block
      Ioss::ElementBlock* eb = region->get_element_block("BLOCK_1");
      EXPECT_TRUE(nullptr != eb);

      int64_t localElementCount = eb->entity_count();
      int64_t globalElementCount = stk::get_global_sum(communicator, localElementCount);

      int64_t goldGlobalElementCount = 8 - i;
      EXPECT_EQ(goldGlobalElementCount, globalElementCount);

      // Should be 1 element field on database named "density";
      EXPECT_EQ(1u, eb->field_count(Ioss::Field::TRANSIENT));
      EXPECT_TRUE(eb->field_exists("density"));

      // Verify that the correct data was written for each step per mesh group.
      double time = i;

      double db_time = region->begin_state(1);
      EXPECT_EQ(time, db_time);

      std::vector<double> field_data;
      eb->get_field_data("density", field_data);
      const double expectedValue = 10.0 * time;
      for (size_t elem = 0; elem < field_data.size(); elem++) {
        EXPECT_NEAR(expectedValue, field_data[elem], 1.0e-6);
      }
      region->end_state(1);
    }
  }

  // ============================================================
  // Cleanup
  unlink(mesh_name.c_str());
  unlink(results_name.c_str());
}

}
