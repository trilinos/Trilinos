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

#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <exception>                    // for exception
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <gtest/gtest.h>
#include <string>                       // for string
#include "Ioss_Field.h"                 // for Field, etc
#include "gtest/gtest.h"                // for AssertHelper, EXPECT_ANY_THROW
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class FieldBase; } }

TEST(StkMeshIoBroker, CheckInvalidCallOrdering)
{
    const std::string outputFilename = "invalid_checks.exo";
    MPI_Comm communicator = MPI_COMM_WORLD;

    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    size_t input_index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(input_index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK, "displacement", 1);
    stk::mesh::put_field_on_mesh(field0, stkMeshMetaData.universal_part(), nullptr);
    stkIo.populate_bulk_data();

    {
      size_t results_output_index = stkIo.create_output_mesh(outputFilename, stk::io::WRITE_RESULTS);

      stk::mesh::FieldBase *field0a = stkMeshMetaData.get_field(stk::topology::NODE_RANK, "displacement");
      stkIo.add_field(results_output_index, *field0a);
      stkIo.add_global(results_output_index, "NotTooLate", "scalar", Ioss::Field::DOUBLE);

      // Global variables defined, process_output_request does not output globals...
      EXPECT_ANY_THROW(stkIo.process_output_request(results_output_index, 0.0));

      stkIo.begin_output_step(results_output_index, 1.0);
      stkIo.write_defined_output_fields(results_output_index);
      stkIo.end_output_step(results_output_index);

      // Try to write a global field after the step has been ended.
//      EXPECT_ANY_THROW(stkIo.write_global(results_output_index, "NotTooLate", 1.0));

      // Try to add a field after output has already been done...
      EXPECT_ANY_THROW(stkIo.add_field(results_output_index, *field0a));

      // Try to add a global field after output has already been done...
      EXPECT_ANY_THROW(stkIo.add_global(results_output_index, "TooLate", "scalar", Ioss::Field::DOUBLE));

      // Try to set the subset selector after output mesh has already been written.
      stk::mesh::Selector selector;
      EXPECT_ANY_THROW(stkIo.set_subset_selector(results_output_index, selector));

      // Try to set the use_nodeset_parts_for_node_fieldssubset selector after output mesh has already been written.
      EXPECT_ANY_THROW(stkIo.use_nodeset_for_block_nodes_fields(results_output_index, true));
      EXPECT_ANY_THROW(stkIo.use_nodeset_for_sideset_nodes_fields(results_output_index, true));

      // Try to call write_defined_output_fields without beginning an output step...
      EXPECT_ANY_THROW(stkIo.write_defined_output_fields(results_output_index));

      // Try to call end_output_step before beginning an output step...
      EXPECT_ANY_THROW(stkIo.end_output_step(results_output_index));

      // Try to call begin_output_step() without calling end_output_step().
      stkIo.begin_output_step(results_output_index, 1.0);
//      EXPECT_ANY_THROW(stkIo.begin_output_step(results_output_index, 1.0));
      stkIo.end_output_step(results_output_index);
    }

    unlink(outputFilename.c_str());
}
